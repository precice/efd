#pragma once

#include "Controller.hpp"

#include "Simulation/Configuration.hpp"
#include "Simulation/Grid.hpp"
#include "Simulation/IfsfdCellAccessor.hpp"
#include "Simulation/IfsfdMemory.hpp"
#include "Simulation/SfsfdCellAccessor.hpp"
#include "Simulation/SfsfdMemory.hpp"

#include <precice/SolverInterface.hpp>

#include <Eigen/Dense>
#include <petscksp.h>

#include <Uni/ExecutionControl/exception>
#include <Uni/Helpers/macros>

#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
//

template <typename TVector>
class RbfBasedInterfaceCell : public InterfaceCell<TVector> {
public:
  using Base = InterfaceCell<TVector>;

  enum {
    Dimensions = Base::Dimensions
  };

  using Vector = typename Base::Vector;

  using Scalar = typename Base::Scalar;

  using InternalIds = Eigen::Matrix<int, Dimensions, 1>;

  RbfBasedInterfaceCell(unsigned const& global_index,
                        Vector const&   position)
    : _globalIndex(global_index),
    _position(position),
    _internalIds(InternalIds::Constant(-1)) {}

  unsigned
  globalIndex() const {
    return _globalIndex;
  }

  Vector const&
  data() const {
    return _data;
  }

  Vector&
  data() {
    return _data;
  }

  Scalar const&
  data(unsigned const& dimension) const {
    return _data(dimension);
  }

  Scalar&
  data(unsigned const& dimension) {
    return _data(dimension);
  }

  bool
  isInternal(unsigned const& dimension) const {
    return _internalIds(dimension) >= 0;
  }

  int const&
  internalIds(unsigned const& dimension) const {
    return _internalIds(dimension);
  }

  int&
  internalIds(unsigned const& dimension) {
    return _internalIds(dimension);
  }

private:
  unsigned _globalIndex;
  Uni_PrivateProperty(Vector,      position);
  Uni_PrivateProperty(InternalIds, internalIds);
  Vector _data;
};

template <typename TVector>
using RbfBasedInterfaceCellSet
        = boost::multi_index_container
          < std::shared_ptr < RbfBasedInterfaceCell < TVector >>,
      boost::multi_index::indexed_by
      < boost::multi_index::ordered_unique
      < boost::multi_index::const_mem_fun
      < RbfBasedInterfaceCell<TVector>,
      unsigned,
      &RbfBasedInterfaceCell<TVector>::globalIndex >> >>;

template <typename TVector>
class RbfBasedIteratorBackEnd :
  public IteratorBackEnd < InterfaceCell < TVector >> {
public:
  using Value = InterfaceCell<TVector>;

  using Base = IteratorBackEnd<Value>;

  using Set = RbfBasedInterfaceCellSet<TVector>;

  using SetIterator
          = typename Set::template nth_index<0>::type::iterator;

  RbfBasedIteratorBackEnd(SetIterator const& it) : _it(it) {}

  RbfBasedIteratorBackEnd(RbfBasedIteratorBackEnd const&) = delete;

  ~RbfBasedIteratorBackEnd() {}

  RbfBasedIteratorBackEnd&
  operator=(RbfBasedIteratorBackEnd const&) = delete;

  std::unique_ptr<Base>
  clone() const {
    return std::unique_ptr<Base>(new RbfBasedIteratorBackEnd(_it));
  }

  bool
  equals(std::unique_ptr<Base> const& other) {
    return _it
           == reinterpret_cast<RbfBasedIteratorBackEnd const*>
           (other.get())->_it;
  }

  Value&
  dereference() const {
    return *_it->get();
  }

  void
  increment() {
    ++_it;
  }

  void
  decrement() {
    --_it;
  }

private:
  SetIterator _it;
};

template <typename TVector>
class RbfBasedIterableBackEnd :
  public IterableBackEnd < InterfaceCell < TVector >> {
public:
  using Vector = TVector;

  using Value = InterfaceCell<Vector>;

  using Base = IterableBackEnd<Value>;

  using BaseIteratorType = typename Base::IteratorType;

  using BaseUniqueIteratorBackEndType
          = typename BaseIteratorType::UniqueBackEndType;

  using IteratorBackEndType = RbfBasedIteratorBackEnd<Vector>;

  using Set = RbfBasedInterfaceCellSet<TVector>;

  RbfBasedIterableBackEnd(Set* set) : _set(set) {}

  RbfBasedIterableBackEnd(RbfBasedIterableBackEnd const&) = delete;

  ~RbfBasedIterableBackEnd() {}

  RbfBasedIterableBackEnd&
  operator=(RbfBasedIterableBackEnd const&) = delete;

  unsigned
  size() {
    return static_cast<unsigned>(_set->size());
  }

  BaseIteratorType
  begin() {
    return BaseIteratorType(
      BaseUniqueIteratorBackEndType(
        new IteratorBackEndType(_set->template get<0>().begin())));
  }

  BaseIteratorType
  end() {
    return BaseIteratorType(
      BaseUniqueIteratorBackEndType(
        new IteratorBackEndType(_set->template get<0>().end())));
  }

  BaseIteratorType
  find(unsigned const& global_index) {
    return BaseIteratorType(
      BaseUniqueIteratorBackEndType(
        new IteratorBackEndType(
          _set->template get<0>().find(global_index))));
  }

private:
  Set* _set;
};

template <typename TSolverTraits>
class RbfBasedController :
  public Controller<typename TSolverTraits::VectorDsType> {
public:
  class LagrangianNode;

public:
  using SolverTraitsType = TSolverTraits;

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

  using BaseType = Controller<VectorDsType>;

  using IterableType = typename BaseType::IterableType;

private:
  using BaseInterfaceCellType = typename BaseType::InterfaceCellType;

  using InterfaceCellType = RbfBasedInterfaceCell<VectorDsType>;

  using BaseIterableBackEndType
          = typename IterableType::BackEndType;

  using BaseSharedIterableBackEndType
          = typename IterableType::SharedBackEndType;

  using IterableBackEndType = RbfBasedIterableBackEnd<VectorDsType>;

  using Set = RbfBasedInterfaceCellSet<VectorDsType>;

public:
  RbfBasedController(
    Configuration const* configuration,
    MemoryType const*    memory) :
    _memory(memory),
    _doComputeSupportRadius(true),
    _doComputeImqShape(true) {
    if (!configuration->isOfType<std::string>(
          "/Ib/Schemes/DirectForcing/RbfBased/SupportRadius")) {
      _doComputeSupportRadius = false;
      _suppportRadius
        = static_cast<ScalarType>(
        configuration->get<long double>(
          "/Ib/Schemes/DirectForcing/RbfBased/SupportRadius"));
    }

    if (!configuration->isOfType<std::string>(
          "/Ib/Schemes/DirectForcing/RbfBased/ImqShape")) {
      _doComputeImqShape = false;
      _imqShape
        = static_cast<ScalarType>(
        configuration->get<long double>(
          "/Ib/Schemes/DirectForcing/RbfBased/ImqShape"));
    }
  }

  void
  initialize(precice::SolverInterface* precice_interface) {
    _preciceInterface = precice_interface;

    _bodyMeshSet.insert(_preciceInterface->getMeshID("BodyMesh"));

    if (_doComputeSupportRadius) {
      _suppportRadius
        = _memory->gridGeometry()->maxCellWidth().maxCoeff();
    }

    if (_doComputeImqShape) {
      _imqShape
        = _memory->parameters()->re()
          / _memory->gridGeometry()->maxCellWidth().maxCoeff();
    }

    logInfo("RBF shape = {1}",      _imqShape);
    logInfo("Support radius = {1}", _suppportRadius);
  }

  void
  precompute() {
    logInfo("Locate Interface Cells has been starting ...");

    auto     mesh_handle = _preciceInterface->getMeshHandle("BodyMesh");
    unsigned body_id     = 0;

    double                               max_temp = -1;
    double                               mix_temp = 100000000000;
    Eigen::Matrix<double, Dimensions, 1> temp_coords2;

    _lagrangians.clear();
    _lagrangians.reserve(mesh_handle.vertices().size());

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

      VectorDsType body_position = temp_coords.template cast<ScalarType>();

      _lagrangians.emplace_back(body_id, body_position);
      ++body_id;
    }

    _velocitySet.clear();
    _forceSet.clear();

    for (unsigned d = 0; d < Dimensions; ++d) {
      _locateCells(d);
    }

    logInfo("Locate Interface Cells has finished");
  }

  void
  processVelocities() {}

  void
  processForces() {
    for (unsigned d = 0; d < Dimensions; ++d) {
      _resolveForces(d);
    }
  }

  IterableType
  getVelocityIterable() {
    using CompoundIterableType
            = CompoundIterableBackEnd<BaseInterfaceCellType>;

    return IterableType(
      BaseSharedIterableBackEndType(
        new CompoundIterableType(
          typename CompoundIterableType::SharedIterableBackEndType(
            new IterableBackEndType(&_velocitySet)),
          typename CompoundIterableType::SharedIterableBackEndType(
            new IterableBackEndType(&_forceSet)))));
  }

  IterableType
  getForceIterable() {
    return IterableType(
      BaseSharedIterableBackEndType(
        new IterableBackEndType(&_forceSet)));
  }

private:
  void
  _locateCells(unsigned const& dimension) {
    _lagrangiansSupports[dimension].clear();

    auto         width = _memory->gridGeometry()->minCellWidth();
    VectorDsType grid_shift;
    grid_shift             = 0.5 * width;
    grid_shift(dimension) += 0.5 * width(dimension);

    VectorDiType support_size
      = (VectorDsType::Constant(_suppportRadius).cwiseQuotient(width)
         + VectorDsType::Ones()).template cast<int>();

    unsigned internal_id = 0;

    for (auto& lagrangian_node : _lagrangians) {
      _lagrangiansSupports[dimension].emplace_back(&lagrangian_node);

      VectorDiType size
        = (lagrangian_node.position() - grid_shift)
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

        auto dist = (lagrangian_node.position() - position).norm();

        if (dist > _suppportRadius) {
          continue;
        }

        unsigned global_index = computeGlobalIndex(global_space_index);

        auto fluid_cell = _insertFluidCell(dimension,
                                           global_index,
                                           internal_id,
                                           position);

        _lagrangiansSupports[dimension].back().addInterfaceCell(fluid_cell);
      }
    }
    _interfaceCellsSize(dimension) = internal_id;
  }

  InterfaceCellType*
  _insertFluidCell(unsigned const&     dimension,
                   unsigned const&     global_index,
                   unsigned&           internal_id,
                   VectorDsType const& position) {
    bool isInternal = (_preciceInterface->inquirePosition(
                         position.template cast<double>().data(),
                         _bodyMeshSet)
                       != precice::constants::positionOutsideOfGeometry());
    auto find_it = _velocitySet.template get<0>().find(global_index);

    if (find_it == _velocitySet.template get<0>().end()) {
      find_it = _forceSet.template get<0>().find(global_index);

      if (find_it == _forceSet.template get<0>().end()) {
        if (isInternal) {
          auto status = _forceSet.emplace(
            std::make_shared<InterfaceCellType>(global_index, position));
          (*status.first)->internalIds(dimension) = internal_id;
          ++internal_id;

          return status.first->get();
        } else {
          auto status = _velocitySet.emplace(
            std::make_shared<InterfaceCellType>(global_index, position));

          return status.first->get();
        }
      } else {
        if (isInternal) {
          if ((*find_it)->internalIds(dimension) == -1) {
            (*find_it)->internalIds(dimension) = internal_id;
            ++internal_id;
          }
        }

        return find_it->get();
      }
    } else {
      if (isInternal) {
        auto status = _forceSet.emplace(*find_it);
        _velocitySet.erase(find_it);
        (*status.first)->internalIds(dimension) = internal_id;
        ++internal_id;

        return status.first->get();
      } else {
        return find_it->get();
      }
    }
  }

  void
  _resolveForces(unsigned const& dimension) {
    using Matrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

    unsigned max_size = 0;
    unsigned mix_size = 1000000000;

    for (auto& supports : _lagrangiansSupports[dimension]) {
      Matrix m(supports.size() + 1, supports.size() + 1);

      max_size
        = std::max(max_size, static_cast<unsigned>(supports.size()));
      mix_size
        = std::min(mix_size, static_cast<unsigned>(supports.size()));

      Vector b(supports.size() + 1);

      for (unsigned j = 0; j < supports.size(); ++j) {
        auto const& fluid_cell = supports[j];

        for (unsigned k = 0; k < j; ++k) {
          auto const& fluid_cell2 = supports[k];
          m(j, k) = computeImqRbf(fluid_cell.cell()->position(),
                                  fluid_cell2.cell()->position());
        }

        for (unsigned k = j + 1; k < supports.size(); ++k) {
          auto const& fluid_cell2 = supports[k];
          m(j, k) = computeImqRbf(fluid_cell.cell()->position(),
                                  fluid_cell2.cell()->position());
        }

        m(j, j)               = computeImqRbf(0.0);
        m(j, supports.size()) = 1.0;
        b(j)                  = computeImqRbf(
          supports.lagrangianNode()->position(),
          fluid_cell.cell()->position());
      }
      m.row(supports.size()).fill(1.0);
      m(supports.size(), supports.size()) = 0.0;
      b(supports.size())                  = 1.0;

      Eigen::ColPivHouseholderQR<Matrix> dec(m);
      Vector                             x = dec.solve(b);

      for (unsigned j = 0; j < supports.size(); ++j) {
        supports[j].weight() = x(j);
        // if (x(j) == 0.0) {
        // }
      }
    }

    logInfo("Max in {1} is {2}", dimension, max_size);
    logInfo("Mix in {1} is {2}", dimension, mix_size);

    /*
     * Solve least-square problem
     */

    // Create right-hand side vector and vector of unknowns
    logInfo("Create b = {1}", _lagrangians.size());
    Vec b, x;
    Vec bt;
    VecCreate(_memory->parallelDistribution()->mpiCommunicator,
              &b);
    VecSetSizes(b, _lagrangians.size(), PETSC_DECIDE);
    VecSetType(b, VECMPI);
    VecCreate(_memory->parallelDistribution()->mpiCommunicator,
              &bt);
    VecSetSizes(bt, _interfaceCellsSize(dimension), PETSC_DECIDE);
    VecSetType(bt, VECMPI);
    logInfo("Create x = {1}", _interfaceCellsSize(dimension));
    VecCreate(_memory->parallelDistribution()->mpiCommunicator,
              &x);
    VecSetSizes(x, _interfaceCellsSize(dimension), PETSC_DECIDE);
    VecSetType(x, VECMPI);

    // Fill right-hand side;

    double temp_row         = 0;
    double temp_norm_base   = 0;
    double temp_norm_base_b = 0;

    for (auto& supports : _lagrangiansSupports[dimension]) {
      PetscScalar value = 0;

      for (auto& fluid_cell : supports) {
        value -= fluid_cell.weight() * fluid_cell.cell()->data(
          dimension);
        temp_norm_base += fluid_cell.cell()->data(dimension)
                          / _memory->timeStepSize();

        temp_norm_base_b += fluid_cell.weight() * fluid_cell.cell()->data(
          dimension)
                            / _memory->timeStepSize();

        if (supports.lagrangianNode()->id() == 1) {
          temp_row += fluid_cell.weight();
        }
      }
      value /= _memory->timeStepSize();

      PetscInt global_index = supports.lagrangianNode()->id();
      VecSetValues(b, 1, &global_index, &value, INSERT_VALUES);
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    logInfo("========= ROW {1}", temp_row);

    // Create matrix
    logInfo("Create matrix");
    Mat global_matrix;
    Mat preconditioning_global_matrix;
    MatCreate(_memory->parallelDistribution()->mpiCommunicator,
              &global_matrix);
    MatSetSizes(global_matrix,
                PETSC_DECIDE,
                PETSC_DECIDE,
                _lagrangiansSupports[dimension].size(),
                _interfaceCellsSize(dimension));
    // MatSetType(global_matrix, MATMPIAIJ);
    MatSetType(global_matrix, MATAIJ);
    MatSetUp(global_matrix);

    // Fill matrix
    // PetscInt end;
    // MatGetOwnershipRange(global_matrix, &first, &end);
    logInfo("Fill matrix");

    double temp_column = 0;

    for (auto const& supports : _lagrangiansSupports[dimension]) {
      PetscInt global_row = supports.lagrangianNode()->id();

      for (auto const& fluid_cell : supports) {
        PetscInt global_column = fluid_cell.cell()->internalIds(dimension);

        if (!fluid_cell.cell()->isInternal(dimension)) {
          continue;
        }
        PetscScalar value = fluid_cell.weight();
        MatSetValues(global_matrix,
                     1, &global_row,
                     1, &global_column,
                     &value, INSERT_VALUES);

        if (global_column == 1) {
          temp_column += value;
        }
      }
    }
    MatAssemblyBegin(global_matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(global_matrix, MAT_FINAL_ASSEMBLY);

    logInfo("========= COLUMN {1}", temp_column);

    KSP solver;
    KSPCreate(_memory->parallelDistribution()->mpiCommunicator,
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
#if ((PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 5))
    KSPSetOperators(solver,
                    // global_matrix,
                    preconditioning_global_matrix,
                    preconditioning_global_matrix
                    // 0
                    );
#else
    KSPSetOperators(solver,
                    // global_matrix,
                    preconditioning_global_matrix,
                    preconditioning_global_matrix,
                    // 0,
                    SAME_PRECONDITIONER);
#endif

    logInfo("KSP Setup");
    KSPSetUp(solver);
    logInfo("KSP Solve");
    MatMultTranspose(global_matrix, b, bt);
    KSPSolve(solver, bt, x);

    MatMult(global_matrix, x, b);
    PetscScalar temp_norm;
    VecNorm(x, NORM_1, &temp_norm);
    logInfo("NOOOOOOOOOOOOORM1 {1}", temp_norm + temp_norm_base);
    VecNorm(b, NORM_1, &temp_norm);
    logInfo("NOOOOOOOOOOOOORM2 {1}", temp_norm + temp_norm_base_b);

    PetscScalar global_norm;
    PetscInt    iteration_number;
    VecNorm(x, NORM_2, &global_norm);
    KSPGetIterationNumber(solver, &iteration_number);

    PetscPrintf(_memory->parallelDistribution()->mpiCommunicator,
                "Norm of vector %G iterations %D\n",
                global_norm,
                iteration_number);

    // Retrieve data
    logInfo("Retrieve data");
    PetscScalar* unknowns;
    VecGetArray(x, &unknowns);

    for (auto& cell : _forceSet) {
      if (cell->isInternal(dimension)) {
        cell->data(dimension) = unknowns[cell->internalIds(dimension)];
      } else {
        cell->data(dimension) = 0.0;
      }
    }

    VecRestoreArray(x, &unknowns);

    logInfo("Clean up");
    VecDestroy(&b);
    VecDestroy(&bt);
    VecDestroy(&x);
    MatDestroy(&global_matrix);
    MatDestroy(&preconditioning_global_matrix);
    KSPDestroy(&solver);
  }

  //

public:
  class LagrangianNode {
public:
    LagrangianNode(unsigned const&     id,
                   VectorDsType const& position)
      : _id(id),
      _position(position) {}

    Uni_PublicProperty(unsigned,     id);
    Uni_PublicProperty(VectorDsType, position);
  };

  class LagrangianNodeSupports {
public:
    class Support {
public:
      Support(InterfaceCellType const* cell,
              ScalarType const&        weight)
        : _cell(cell),
        _weight(weight) {}

      Uni_PublicProperty(InterfaceCellType const*, cell);
      Uni_PublicProperty(ScalarType,               weight);
    };
    using Supports = std::vector<Support>;

    LagrangianNodeSupports(LagrangianNode* lagrangian_node)
      : _lagrangianNode(lagrangian_node) {}

    typename Supports::const_iterator
    begin() const {
      return _supports.begin();
    }

    typename Supports::iterator
    begin() {
      return _supports.begin();
    }

    typename Supports::const_iterator
    end() const {
      return _supports.end();
    }

    typename Supports::iterator
    end() {
      return _supports.end();
    }

    std::size_t
    size() const {
      return _supports.size();
    }

    void
    addInterfaceCell(InterfaceCellType* cell) {
      _supports.emplace_back(cell, 0.0);
    }

    Support&
    operator[](std::size_t const& index) {
      return _supports[index];
    }

    Uni_PublicProperty(LagrangianNode*, lagrangianNode);
    Uni_PublicProperty(Supports,        supports);
  };

  class AccessorTraits {
public:
    enum {
      Dimensions = RbfBasedController::Dimensions
    };

    using Type
            = Uni::StructuredGrid::Basic::GlobalMultiIndex<AccessorTraits>;

    using GridType
            = Uni::StructuredGrid::Basic::Grid<Type>;
  };

private:
  ScalarType
  computeImqRbf(VectorDsType one, VectorDsType two) const {
    auto radius = (one - two).norm();

    return computeImqRbf(radius);
  }

  ScalarType
  computeImqRbf(ScalarType radius) const {
    return 1.0 / std::sqrt(1.0 + _imqShape * _imqShape * radius * radius);
  }

  unsigned
  computeGlobalIndex(VectorDiType const& global_space_index) const {
    auto temp = global_space_index
                + _memory->grid()->innerGrid.leftIndent();

    unsigned global_index = temp(Dimensions - 1);

    for (int d = Dimensions - 2; d >= 0; --d) {
      global_index *= _memory->grid()->innerGrid.size(d);
      global_index += temp(d);
    }

    return global_index;
  }

  precice::SolverInterface* _preciceInterface;
  MemoryType const*         _memory;

  std::set<int> _bodyMeshSet;

  bool       _doComputeSupportRadius;
  ScalarType _suppportRadius;
  bool       _doComputeImqShape;
  ScalarType _imqShape;

  Eigen::Matrix<unsigned, Dimensions, 1> _interfaceCellsSize;

  Set _velocitySet;
  Set _forceSet;

  std::vector<LagrangianNode> _lagrangians;
  std::array<std::vector<LagrangianNodeSupports>,
             Dimensions> _lagrangiansSupports;
};
}
}
}
