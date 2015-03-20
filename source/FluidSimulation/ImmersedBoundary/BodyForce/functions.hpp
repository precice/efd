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
template <typename TCellAccessor>
bool
computeCellForce(TCellAccessor const&                      accessor,
                 typename TCellAccessor::ScalarType const& re,
                 typename TCellAccessor::VectorDsType&     force) {
  // typedef Eigen::Matrix<typename TCellAccessor::ScalarType,
  //                       TCellAccessor::Dimensions,
  //                       TCellAccessor::Dimensions> Matrix;

  // typedef typename TCellAccessor::VectorDsType Vector;

  // bool isCurrentOutside  = true;
  // bool isNeighborsInside = false;

  // for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
  //   if (accessor.positionInRespectToGeometrys(d)
  //       != precice::constants::positionOutsideOfGeometry()) {
  //     isCurrentOutside = false;
  //     break;
  //   }

  //   if (isNeighborsInside) {
  //     continue;
  //   }

  //   for (int d2 = 0; d2 < TCellAccessor::Dimensions; ++d2) {
  //     for (int d3 = 0; d3 < 2; ++d3) {
  //       if (accessor.relativeCell(d2, d3)->positionInRespectToGeometrys(d)
  //           != precice::constants::positionOutsideOfGeometry()) {
  //         isNeighborsInside = true;
  //         break;
  //       }
  //     }
  //   }
  // }

  // if (!(isCurrentOutside && isNeighborsInside)) {
  //   return false;
  // }

  // Matrix matrix;

  // for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
  //   for (int d2 = 0; d2 < TCellAccessor::Dimensions; ++d2) {
  //     matrix(d, d2) = dudx(
  //       accessor.relativeCell(d2, 0)->velocity(d),
  //       accessor.relativeCell(d2, 1)->velocity(d),
  //       accessor.width(d2, 0) (d),
  //       accessor.width(d2, 1) (d));
  //   }
  // }

  // matrix = 1.0 / re * (matrix + matrix.transpose());

  // matrix -= accessor.pressure() * Matrix::Identity();

  // for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
  //   for (int d2 = 0; d2 < 2; ++d2) {
  //     for (int d3 = 0; d3 < TCellAccessor::Dimensions; ++d3) {
  //       if (accessor.relativeCell(d, d2)->positionInRespectToGeometrys(d3)
  //           != precice::constants::positionOutsideOfGeometry()) {
  //         Vector normal = Vector::Zero();

  //         if (d2 == 0) {
  //           normal(d) = 1.0;
  //         } else {
  //           normal(d) = -1.0;
  //         }

  //         for (int d4 = 0; d4 < TCellAccessor::Dimensions; ++d4) {
  //           if (d4 != d) {
  //             normal(d) *= accessor.width() (d4);
  //           }
  //         }
  //         force += matrix * normal;
  //         // logInfo("Got it {1} {2} {3}", matrix, normal, force);
  //         break;
  //       }
  //     }
  //   }
  // }

  //
  return true;
}
}
}
}
}

#endif
