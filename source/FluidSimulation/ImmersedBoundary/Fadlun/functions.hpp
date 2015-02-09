#ifndef FsiSimulation_FluidSimulation_ImmersedBoundary_Fadlun_functions_hpp
#define FsiSimulation_FluidSimulation_ImmersedBoundary_Fadlun_functions_hpp

#include <precice/SolverInterface.hpp>

#include <Uni/Logging/macros>

#include <Eigen/Core>

#include <array>
#include <cmath>
#include <limits>

namespace FsiSimulation {
namespace FluidSimulation {
namespace ImmersedBoundary {
namespace Fadlun {
template <typename TCellAccessor,
          int TDimensions>
void
computeDistances(TCellAccessor const&      accessor,
                 precice::SolverInterface* preciceInterface) {
  typedef Eigen::Matrix<double, TDimensions, 1> Position;
  typedef std::array<Position, TDimensions>     VelocityPositions;

  auto currentPosition = accessor.currentPosition();
  auto currentWidth    = accessor.currentWidth();

  VelocityPositions positions;

  for (int d = 0; d < TDimensions; ++d) {
    positions[d]     = currentPosition.template cast<double>();
    positions[d]    += 0.5 * currentWidth.template cast<double>();
    positions[d](d) += 0.5 * (double)currentWidth(d);
    // logInfo("1 {1}", positions[d].transpose());
  }

  for (int d = 0; d < TDimensions; ++d) {
    auto position = preciceInterface->inquirePosition(
      positions[d].data(),
      std::set<int>({ preciceInterface->getMeshID("Body") }));
    accessor.currentCell()->positions(d) = position;
  }

  auto position = (accessor.currentPosition() + 0.5 *
                   accessor.currentWidth()).eval();
  accessor.currentCell()->position()
    = preciceInterface->inquirePosition(
    position.data(),
    std::set<int>({ preciceInterface->getMeshID("Body") }));
}

template <typename TCellAccessor,
          typename TScalar,
          int TDimensions>
inline void
computeVelocity(
  TCellAccessor const&                        accessor,
  typename TCellAccessor::VectorDsType const& boundaryVelocity,
  int const&                                  d,
  int const&                                  d2) {
  TScalar const width =
    0.5 * accessor.currentWidth() (d)
    + 0.5 * accessor.relativeWidth(d, d2) (d);
  TScalar const b1 = accessor.currentCell()->distances(d);
  TScalar const b2 = accessor.relativeCell(d, d2)->distances(d);
  TScalar const a1 = width * b1 / (b2 + b1);
  TScalar const a2 = width - a1;

  if (a1 < 0 || a2 < 0) {
    logInfo("Less than 0");

    return;
  }
  TScalar const width2 =
    0.5 * accessor.currentWidth() (d)
    + 0.5 * accessor.relativeWidth(d, !d2) (d);
  accessor.currentCell()->velocity(d)
  // = boundaryVelocity(d);
    = (a1 * boundaryVelocity(d)
       + width2 * accessor.relativeCell(d, !d2)->velocity(d))
      / (a1 + width2);
}

template <typename TCellAccessor,
          typename TScalar,
          int TDimensions>
bool
treatBoundary(TCellAccessor const&      accessor,
              precice::SolverInterface* preciceInterface,
              TScalar const&            dt,
              int const&                d) {
  typedef typename TCellAccessor::CellType::Velocity Velocity;
  typedef typename TCellAccessor::VectorDsType       VectorDs;

  Velocity const boundaryVelocity = Velocity::Zero();

  if (accessor.currentCell()->positions(d)
      == precice::constants::positionOutsideOfGeometry()) {
    for (int d2 = 0; d2 < 2; ++d2) {
      if (accessor.relativeCell(d, d2)->positions(d)
          == precice::constants::positionInsideOfGeometry()) {
        if (accessor.relativeCell(d, !d2)->positions(d)
            == precice::constants::positionOutsideOfGeometry()) {
          computeVelocity<TCellAccessor, TScalar, TDimensions>(
            accessor,
            boundaryVelocity,
            d, d2);

          return true;
        } else {
          logError("This is bad {1}", accessor.indexValues().transpose());
          accessor.currentCell()->velocity(d)
            = boundaryVelocity(d);

          return true;
        }
      }
    }
  } else if (accessor.currentCell()->positions(d)
             == precice::constants::positionInsideOfGeometry()) {
    accessor.currentCell()->velocity(d) = boundaryVelocity(d);

    return true;

    for (int d2 = 0; d2 < 2; ++d2) {
      if (accessor.relativeCell(d, d2)->positions(d)
          == precice::constants::positionOutsideOfGeometry()) {
        if (accessor.relativeCell(d, !d2)->positions(d)
            == precice::constants::positionInsideOfGeometry()) {
          computeVelocity<TCellAccessor, TScalar, TDimensions>(
            accessor,
            boundaryVelocity,
            d, d2);

          return true;
        } else {
          logWarning("Gottya {1}", accessor.indexValues().transpose());
          accessor.currentCell()->velocity(d)
            = boundaryVelocity(d);

          return true;
        }
      }
    }
  } else if (accessor.currentCell()->positions(d)
             == precice::constants::positionOnGeometry()) {
    accessor.currentCell()->velocity(d) = boundaryVelocity(d);

    return true;
  }

  return false;
}
}
}
}
}

#endif
