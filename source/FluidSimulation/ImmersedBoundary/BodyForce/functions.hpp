#ifndef FsiSimulation_FluidSimulation_ImmersedBoundary_BodyForce_functions_hpp
#define FsiSimulation_FluidSimulation_ImmersedBoundary_BodyForce_functions_hpp

#include "FluidSimulation/functions.hpp"

#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <precice/SolverInterface.hpp>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
namespace BodyForce {
template <typename TCellAccessor,
          int TDimensions>
bool
computeCellForce(TCellAccessor const&                        accessor,
                 typename TCellAccessor::CellType::Scalar const& re,
                 typename TCellAccessor::CellType::VectorDs&     force) {
  typedef Eigen::Matrix<typename TCellAccessor::CellType::Scalar,
                        TDimensions, TDimensions> Matrix;
  typedef typename TCellAccessor::CellType::VectorDs Vector;

  bool isCurrentOutside  = true;
  bool isNeighborsInside = false;

  for (int d = 0; d < TDimensions; ++d) {
    if (accessor.currentCell()->positions(d)
        != precice::constants::positionOutsideOfGeometry()) {
      isCurrentOutside = false;
      break;
    }

    if (isNeighborsInside) {
      continue;
    }

    for (int d2 = 0; d2 < TDimensions; ++d2) {
      for (int d3 = 0; d3 < 2; ++d3) {
        if (accessor.relativeCell(d2, d3)->positions(d)
            != precice::constants::positionOutsideOfGeometry()) {
          isNeighborsInside = true;
          break;
        }
      }
    }
  }

  if (!(isCurrentOutside && isNeighborsInside)) {
    return false;
  }

  Matrix matrix;

  for (int d = 0; d < TDimensions; ++d) {
    for (int d2 = 0; d2 < TDimensions; ++d2) {
      matrix (d, d2) = dudx(
        accessor.relativeCell(d2, 0)->velocity(d),
        accessor.relativeCell(d2, 1)->velocity(d),
        accessor.relativeWidth(d2, 0) (d),
        accessor.relativeWidth(d2, 1) (d));
    }
  }

  matrix = 1.0 / re * (matrix + matrix.transpose());

  matrix -= accessor.currentCell()->pressure() * Matrix::Identity();

  for (int d = 0; d < TDimensions; ++d) {
    for (int d2 = 0; d2 < 2; ++d2) {
      for (int d3 = 0; d3 < TDimensions; ++d3) {
        if (accessor.relativeCell(d, d2)->positions(d3)
            != precice::constants::positionOutsideOfGeometry()) {
          Vector normal = Vector::Zero();

          if (d2 == 0) {
              normal(d) = 1.0;
          } else {
              normal(d) = -1.0;
          }

          for (int d4 = 0; d4 < TDimensions; ++d4) {
            if (d4 != d) {
              normal(d) *= accessor.currentWidth() (d4);
            }
          }
          force += matrix * normal;
          //logInfo("Got it {1} {2} {3}", matrix, normal, force);
          break;
        }
      }
    }
  }

  //
  return true;
}
}
}
}
}

#endif
