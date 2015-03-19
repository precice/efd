#ifndef \
  FsiSimulation_FluidSimulation_ImmersedBoundary_FeedbackForcing_functions_hpp
#define \
  FsiSimulation_FluidSimulation_ImmersedBoundary_FeedbackForcing_functions_hpp

#include <precice/SolverInterface.hpp>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
namespace FeedbackForcing {
template <typename TCellAccessor,
          int TDimensions>
bool
treatBoundary(TCellAccessor const&                            accessor,
              typename TCellAccessor::CellType::Scalar const& alpha,
              int const&                                      d) {
  typedef typename TCellAccessor::CellType::Velocity Velocity;
  typedef typename TCellAccessor::CellType::VectorDs VectorDs;

  Velocity const boundaryVelocity = Velocity::Zero();

  if (accessor.positionInRespectToGeometrys(d)
      != precice::constants::positionOutsideOfGeometry()) {
    for (int d2 = 0; d2 < 2; ++d2) {
      if (accessor.relativeCell(d, d2)->positionInRespectToGeometrys(d)
          == precice::constants::positionOutsideOfGeometry()) {
        accessor.fgh(d) += alpha
                                          * accessor.velocity(d);

        return true;
      }
    }
  }

  return false;
}
}
}
}
}

#endif
