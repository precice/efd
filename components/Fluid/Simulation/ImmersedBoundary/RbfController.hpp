#include "Controller.hpp"

#include <Eigen/Dense>

#include <petscksp.h>

#include <list>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
namespace Private {
template <typename TVector>
struct RbfBodyInterfaceCell : public  BodyInterfaceCell<TVector> {
  using Vector = TVector;

  using BaseType = BodyInterfaceCell<TVector>;

  enum {
    Dimensions = Vector::RowsAtCompileTime
  };

  RbfBodyInterfaceCell() : BaseType(0) {}

  Vector velocity;
};
}
template <typename TSolverTraits>
class RbfController : public Controller<TSolverTraits> {
public:
  using SolverTraitsType = TSolverTraits;

  using BaseType = Controller<SolverTraitsType>;

  using BodyInterfaceCellType = typename BaseType::BodyInterfaceCellType;

  using VertexIdMap = typename BaseType::VertexIdMap;

  enum {
    Dimensions = SolverTraitsType::Dimensions
  };

  using GridGeometryType = typename SolverTraitsType::GridGeometryType;

  using ParametersType = typename SolverTraitsType::ParametersType;

  using ParallelDistributionType
          = typename SolverTraitsType::ParallelDistributionType;

  using MemoryType = typename SolverTraitsType::MemoryType;

  using CellAccessorType = typename SolverTraitsType::CellAccessorType;

  using GridType = typename SolverTraitsType::GridType;

  using VectorDsType = typename SolverTraitsType::VectorDsType;

  using VectorDiType = typename SolverTraitsType::VectorDiType;

  using ScalarType = typename SolverTraitsType::ScalarType;

private:
  class BodyCell;
  class FluidCell;

public:
  RbfController() : BaseType() {}

  void
  initialize(precice::SolverInterface* preciceInterface,
             MemoryType const*         memory) {
    this->BaseType::initialize(preciceInterface, memory);
    _rbfShape
      = this->_memory->parameters()->re()
        / this->_memory->gridGeometry()->maxCellWidth().maxCoeff();
    _suppportRadius
      = this->_memory->gridGeometry()->maxCellWidth().maxCoeff();
    logInfo("RBF shape = {1}",      _rbfShape);
    logInfo("Support radius = {1}", _suppportRadius);
  }

  void
  computePositionInRespectToGeometry(CellAccessorType const& accessor) const {
    this->BaseType::computePositionInRespectToGeometry(accessor);
  }

  void
  createFluidMeshVertex(CellAccessorType const& accessor) {
    set_cell_neighbors_along_geometry_interface(
      accessor,
      this->_preciceInterface,
      std::set<int>({ this->_bodyMeshId }),
      5);

    bool doAdd = validate_layer_number(accessor,
                                       5,
                                       5);

    if (doAdd) {
      // logInfo("{1} | d = {2}", accessor.index().transpose(), distance);
      // VectorDsType position = accessor.pressurePosition();

      this->_vertexIds[0].insert(
        std::make_pair(
          accessor.globalIndex(),
          std::unique_ptr<BodyInterfaceCellType>(
            new Private::RbfBodyInterfaceCell<VectorDsType>())));
    }
  }

  std::pair<typename VertexIdMap::iterator, bool>
  doesVertexExist(CellAccessorType const& accessor) {
    return this->BaseType::doesVertexExist(accessor);
  }

  void
  setFluidVelocity(typename VertexIdMap::iterator const& it,
                   VectorDsType const&                   velocity) {
    reinterpret_cast<Private::RbfBodyInterfaceCell<VectorDsType>*>
    (it->second.get())->velocity = velocity;
  }

  void
  writeFluidVelocities() {}

  void
  mapData(ScalarType const& time_step_size) {
    ((void)time_step_size);

    for (unsigned d = 0; d < Dimensions; ++d) {
      solveForces(d);
    }
  }

  void
  solveForces(unsigned const& dimension) {
    auto         width = this->_memory->gridGeometry()->minCellWidth();
    VectorDsType grid_shift;
    grid_shift             = 0.5 * width;
    grid_shift(dimension) += 0.5 * width(dimension);

    VectorDiType support_size
      = (VectorDsType::Constant(_suppportRadius).cwiseQuotient(width)
         + VectorDsType::Ones()).template cast<int>();

    auto mesh_handle = this->_preciceInterface->getMeshHandle("BodyMesh");

    unsigned body_id     = 0;
    unsigned internal_id = 0;

    double                               max_temp = -1;
    double                               mix_temp = 100000000000;
    Eigen::Matrix<double, Dimensions, 1> temp_coords2;

    for (auto vertex = mesh_handle.vertices().begin();
         vertex != mesh_handle.vertices().end();
         vertex++) {
      Eigen::Matrix<double, Dimensions, 1> temp_coords(vertex.vertexCoords());

      if (max_temp != -1) {
        max_temp = std::max(max_temp,
                            (temp_coords - temp_coords2).norm());
        mix_temp = std::min(mix_temp,
                            (temp_coords - temp_coords2).norm());
      } else {
        max_temp = 0;
      }

      temp_coords2 = temp_coords;

      // Add new vertex

      VectorDsType body_position = temp_coords.template cast<ScalarType>();

      _cells.emplace_back(body_id, body_position);
      ++body_id;

      VectorDiType size
        = (body_position - grid_shift)
          .cwiseQuotient(width).template cast<int>();

      // Compute search boundaries for the fluid cells
      // around the body vertex

      VectorDiType search_offset
        = size - support_size - 0 * VectorDiType::Ones();
      VectorDiType search_size
        = size + support_size + 2 * VectorDiType::Ones();

      typename AccessorTraits::GridType grid;
      grid.initialize(search_size - search_offset);

      for (auto const& accessor : grid) {
        VectorDiType global_space_index = search_offset + accessor.index();

        VectorDsType position
          = global_space_index.template cast<ScalarType>().cwiseProduct(width)
            + grid_shift;

        auto dist = (body_position - position).norm();

        if (dist > _suppportRadius) {
          continue;
        }

        _list.emplace_back(
          position,
          computeGlobalIndex(global_space_index));

        _cells.back().cells.emplace_back(&_list.back());

        auto it = _map.insert(
          std::make_pair(_list.back().globalIndex, &_list.back()));

        if (!it.second) {
          _list.back().internalId = it.first->second->internalId;
        } else {
          if (this->_preciceInterface->inquirePosition(
                _list.back().position.template cast<double>().data(),
                std::set<int>({ this->_bodyMeshId }))
              != precice::constants::positionOutsideOfGeometry()) {
            _list.back().internalId = internal_id;
            ++internal_id;
          } else {
            _list.back().internalId = -1;
          }
        }
      }
    }
    logInfo("MaX ========= {1}", max_temp);
    logInfo("MaX ========= {1}", mix_temp);

    using Matrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

    unsigned max_size = 0;

    for (auto const& body_cell : _cells) {
      Matrix m(body_cell.cells.size() + 1, body_cell.cells.size() + 1);
      Vector rbf(body_cell.cells.size() + 1);

      max_size
        = std::max(max_size, static_cast<unsigned>(body_cell.cells.size()));

      for (unsigned j = 0; j < body_cell.cells.size(); ++j) {
        auto const& fluid_cell = *body_cell.cells[j];

        for (unsigned k = 0; k < j; ++k) {
          auto const& fluid_cell2 = *body_cell.cells[k];
          m(j, k) = computeImqRbf(fluid_cell.position, fluid_cell2.position);
        }

        for (unsigned k = j + 1; k < body_cell.cells.size(); ++k) {
          auto const& fluid_cell2 = *body_cell.cells[k];
          m(j, k) = computeImqRbf(fluid_cell.position, fluid_cell2.position);
        }

        m(j, j)                      = computeImqRbf(0.0);
        m(j, body_cell.cells.size()) = 1.0;
        rbf(j)                       = computeImqRbf(body_cell.position,
                                                     fluid_cell.position);

        if ((body_cell.position - fluid_cell.position).norm()
            < 0.1 * _suppportRadius) {
          logInfo("Issue | {1} | {2}",
                  (body_cell.position - fluid_cell.position).norm(),
                  computeImqRbf((body_cell.position -
                                 fluid_cell.position).norm()));
        }
      }
      m.row(body_cell.cells.size()).fill(1.0);
      m(body_cell.cells.size(), body_cell.cells.size()) = 0.0;
      rbf(body_cell.cells.size())                       = 1.0;

      Eigen::ColPivHouseholderQR<Matrix> dec(m);

      for (unsigned j = 0; j < body_cell.cells.size(); ++j) {
        Vector b(body_cell.cells.size() + 1);
        b.fill(0.0);
        b(j) = 1.0;
        Vector x = dec.solve(b);

        auto& fluid_cell = body_cell.cells[j];
        fluid_cell->w = x.cwiseProduct(rbf).sum();
      }
    }

    logInfo("Max in {1} is {2}", dimension, max_size);

    /*
     * Solver least-square problem
     */

    // Create right-hand side vector and vector of unknowns
    logInfo("Create b = {1}", body_id);
    Vec b, x;
    Vec bt;
    VecCreate(this->_memory->parallelDistribution()->mpiCommunicator,
              &b);
    VecSetSizes(b, body_id, PETSC_DECIDE);
    VecSetType(b, VECMPI);
    VecCreate(this->_memory->parallelDistribution()->mpiCommunicator,
              &bt);
    VecSetSizes(bt, internal_id, PETSC_DECIDE);
    VecSetType(bt, VECMPI);
    logInfo("Create x = {1}", internal_id);
    VecCreate(this->_memory->parallelDistribution()->mpiCommunicator,
              &x);
    VecSetSizes(x, internal_id, PETSC_DECIDE);
    VecSetType(x, VECMPI);

    // Fill right-hand side;
    // PetscInt first;
    // VecGetOwnershipRange(b, &first, PETSC_NULL);
    // PetscInt localsize;
    // VecGetLocalSize(b, &localsize);
    // logInfo("Fill RHS {1} {2}", first, localsize);

    for (auto const& body_cell : _cells) {
      PetscScalar value = 0;

      for (auto const& fluid_cell : body_cell.cells) {
        value -= fluid_cell->w * getVelocity(dimension, fluid_cell);
      }
      value /= this->_memory->timeStepSize();

      PetscInt global_index = body_cell.id;
      VecSetValues(b, 1, &global_index, &value, INSERT_VALUES);
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    // Create matrix
    logInfo("Create matrix");
    Mat global_matrix;
    Mat preconditioning_global_matrix;
    MatCreate(this->_memory->parallelDistribution()->mpiCommunicator,
              &global_matrix);
    MatSetSizes(global_matrix,
                PETSC_DECIDE,
                PETSC_DECIDE,
                body_id,
                internal_id);
    // MatSetType(global_matrix, MATMPIAIJ);
    MatSetType(global_matrix, MATAIJ);
    MatSetUp(global_matrix);

    // Fill matrix
    // PetscInt end;
    // MatGetOwnershipRange(global_matrix, &first, &end);
    logInfo("Fill matrix");

    for (auto const& body_cell : _cells) {
      PetscInt global_row = body_cell.id;

      for (auto const& fluid_cell : body_cell.cells) {
        if (!fluid_cell->isInternal()) {
          continue;
        }
        PetscInt    global_column = fluid_cell->internalId;
        PetscScalar value         = fluid_cell->w;
        MatSetValues(global_matrix,
                     1, &global_row,
                     1, &global_column,
                     &value, INSERT_VALUES);
      }
    }
    MatAssemblyBegin(global_matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(global_matrix, MAT_FINAL_ASSEMBLY);

    KSP solver;
    KSPCreate(this->_memory->parallelDistribution()->mpiCommunicator,
              &solver);
    KSPSetType(solver, KSPGMRES);
    PC preconditioner;
    KSPGetPC(solver, &preconditioner);
    // PCSetType(preconditioner, PCNONE);
    PCSetType(preconditioner, PCLU);
    KSPSetTolerances(solver,
                     1e-21,
                     1e-25,
                     PETSC_DEFAULT,
                     10000);
    MatTransposeMatMult(global_matrix,
                        global_matrix,
                        MAT_INITIAL_MATRIX,
                        PETSC_DEFAULT,
                        &preconditioning_global_matrix);
    KSPSetOperators(solver,
                    // global_matrix,
                    preconditioning_global_matrix,
                    preconditioning_global_matrix,
                    // 0,
                    SAME_PRECONDITIONER);
    logInfo("KSP Setup");
    KSPSetUp(solver);
    logInfo("KSP Solve");
    MatMultTranspose(global_matrix, b, bt);
    KSPSolve(solver, bt, x);

    PetscScalar global_norm;
    PetscInt    iteration_number;
    VecNorm(x, NORM_2, &global_norm);
    KSPGetIterationNumber(solver, &iteration_number);

    PetscPrintf(this->_memory->parallelDistribution()->mpiCommunicator,
                "Norm of vector %G iterations %D\n",
                global_norm,
                iteration_number);

    // Retrieve data
    logInfo("Retrieve data");
    PetscScalar* unknowns;
    VecGetArray(x, &unknowns);

    for (auto const& cell : _map) {
      if (!cell.second->isInternal()) {
        continue;
      }
      setForce(dimension, cell.second, unknowns);
    }

    VecRestoreArray(x, &unknowns);

    logInfo("Clean up");
    VecDestroy(&b);
    VecDestroy(&bt);
    VecDestroy(&x);
    MatDestroy(&global_matrix);
    MatDestroy(&preconditioning_global_matrix);
    KSPDestroy(&solver);
    _cells.clear();
    _map.clear();
    _list.clear();
  }

  void
  readFluidForces() {}

  VectorDsType
  getFluidForce(typename VertexIdMap::iterator const& it) {
    return it->second->data;
  }

private:
  ScalarType
  computeImqRbf(VectorDsType one, VectorDsType two) const {
    auto radius = (one - two).norm();

    return computeImqRbf(radius);
  }

  ScalarType
  computeImqRbf(ScalarType radius) const {
    return 1.0 / std::sqrt(1.0 + _rbfShape * _rbfShape * radius * radius);
  }

  unsigned
  computeGlobalIndex(VectorDiType const& global_space_index) const {
    auto temp = global_space_index
                + this->_memory->grid()->innerGrid.leftIndent();

    unsigned global_index = temp(Dimensions - 1);

    for (int d = Dimensions - 2; d >= 0; --d) {
      global_index *= this->_memory->grid()->innerGrid.size(d);
      global_index += temp(d);
    }

    return global_index;
  }

  ScalarType
  getVelocity(unsigned const& dimension, FluidCell const* cell) const {
    auto find_it = this->_vertexIds[0].find(cell->globalIndex);

    if (find_it == this->_vertexIds[0].end()) {
      throwException("Did not find a fluid cell with index {1}",
                     cell->globalIndex);

      return 0.0;
    }

    return reinterpret_cast<Private::RbfBodyInterfaceCell<VectorDsType>*>
           (find_it->second.get())->velocity(dimension);
  }

  void
  setForce(unsigned const&    dimension,
           FluidCell const*   cell,
           PetscScalar const* unknowns) {
    PetscScalar petsc_force = unknowns[cell->internalId];
    // ((void)petsc_force);

    auto find_it = this->_vertexIds[0].find(cell->globalIndex);

    if (find_it == this->_vertexIds[0].end()) {
      throwException("Did not find a fluid cell with index {1}",
                     cell->globalIndex);
    }
    find_it->second->data(dimension) = petsc_force;
  }

  ScalarType _rbfShape;
  ScalarType _suppportRadius;

  std::vector<BodyCell>          _cells;
  std::list<FluidCell>           _list;
  std::map<unsigned, FluidCell*> _map;

private:
  //

  class FluidCell {
public:
    FluidCell(VectorDsType const& position_,
              unsigned const&     global_index)
      : position(position_),
      globalIndex(global_index) {}

    bool
    isInternal() const {
      return internalId >= 0;
    }

    int          internalId;
    VectorDsType position;
    unsigned     globalIndex;

    int rank;

    ScalarType w;
  };

  class BodyCell {
public:
    BodyCell(unsigned const& id_,
             VectorDsType    position_)
      : id(id_),
      position(position_) {}

    unsigned                id;
    VectorDsType            position;
    std::vector<FluidCell*> cells;
  };

  class AccessorTraits {
public:
    enum {
      Dimensions = RbfController::Dimensions
    };

    using Type
            = Uni::StructuredGrid::Basic::GlobalMultiIndex<AccessorTraits>;

    using GridType
            = Uni::StructuredGrid::Basic::Grid<Type>;
  };
};
}
}
}
