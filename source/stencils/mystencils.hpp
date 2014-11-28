#ifndef FsiSimulation_mystencils_hpp
#define FsiSimulation_mystencils_hpp

#include <Eigen/Core>

#include <mpi.h>
#include <petsc.h>

#include <limits>

namespace FsiSimulation {
template <typename Scalar, int D>
using ScalarVector =  Eigen::Matrix<Scalar, D, 1>;

template <typename Scalar>
inline void
mpiAllReduceMin(Scalar& from, Scalar& to) {}

template <>
inline void
mpiAllReduceMin<float
                >(float& from, float& to) {
  MPI_Allreduce(&from, &to, 1, MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);
}

template <>
inline void
mpiAllReduceMin<double
                >(double& from, double& to) {
  MPI_Allreduce(&from, &to, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
}

template <typename Scalar, int D>
inline Scalar
setTimeStep(Scalar const& re,
            Scalar const& tau,
            ScalarVector<Scalar, D> const& minCellSize,
            ScalarVector<Scalar, D> const& maxVelocity) {
  Scalar localMin, globalMin;
  Scalar factor = 0.0;

  auto square = minCellSize.dot(minCellSize);

  for (auto const& value : square) {
    factor += 1.0 / value;
  }

  auto inverseMaxVelocity(maxVelocity);

  for (auto& value : inverseMaxVelocity) {
    value = 1.0 / value;
  }

  localMin = std::min(re / (2.0 * factor),
                      inverseMaxVelocity.minCoeff());

  globalMin = std::numeric_limits<Scalar>::max();
  mpiAllReduceMin(localMin, globalMin);

  factor  = globalMin;
  factor *= tau;

  return factor;
}

template <typename Scalar>
inline Scalar
dudx(Scalar const& currentU,
     Scalar const& leftU,
     Scalar const& currentX) {
  return (currentU - leftU) / currentX;
}

template <typename Scalar>
inline Scalar
d2udx2(Scalar const& currentU,
       Scalar const& leftU,
       Scalar const& rightU,
       Scalar const& currentX,
       Scalar const& rightX) {
  Scalar const factor = currentX + rightX;

  return 2.0 * (rightU / (rightX * factor) -
                currentU / (currentX * rightX) +
                leftU / (currentX * factor));
}

template <typename Scalar>
inline Scalar
d2udy2(Scalar const& currentU,
       Scalar const& bottomU,
       Scalar const& topU,
       Scalar const& currentY,
       Scalar const& bottomY,
       Scalar const& topY) {
  Scalar const factorBottom = 0.5 * (currentY + bottomY);
  Scalar const factorTop    = 0.5 * (currentY + topY);
  Scalar const factor       = factorBottom + factorTop;

  return 2.0 * (topU / (factorTop * factor) -
                currentU / (factorBottom * factorTop) +
                bottomU / (factorBottom * factor));
}

template <typename Scalar>
inline Scalar
duvdy(Scalar const& currentU,
      Scalar const& currentV,
      Scalar const& bottomU,
      Scalar const& bottomV,
      Scalar const& rightBottomV,
      Scalar const& rightV,
      Scalar const& topU,
      Scalar const& currentY,
      Scalar const& currentX,
      Scalar const& bottomY,
      Scalar const& topY,
      Scalar const& rightX,
      Scalar const& gamma) {
  // distance of corner points in x-direction from center v-value
  Scalar const yCurrent = 0.5 * currentY;
  // distance between center and west v-value
  Scalar const yBottom = 0.5 * (currentY + bottomY);
  // distance between center and east v-value
  Scalar const yTop = 0.5 * (currentY + topY);
  // distance of center u-value from upper edge of cell
  Scalar const xCurrent = 0.5 * currentX;
  // distance of north and center u-value
  Scalar const xRight = 0.5 * (currentX + rightX);

  Scalar const secondOrder = (((xRight - xCurrent) / xRight * currentV +
                               xCurrent / xRight * rightV) *
                              ((yTop - yCurrent) / yTop * currentU +
                               yCurrent / yTop * topU) -
                              ((xRight - xCurrent) / xRight * bottomV +
                               xCurrent / xRight * rightBottomV) *
                              ((yBottom - yCurrent) / yBottom * currentU +
                               yCurrent / yBottom * bottomU)
                              ) / (2.0 * yCurrent);

  Scalar const kr = (xRight - xCurrent) / xRight * currentV +
                    xCurrent / xRight * rightV;
  Scalar const kl = (xRight - xCurrent) / xRight * bottomV +
                    xCurrent / xRight * rightBottomV;

  Scalar const firstOrder = 1.0 / (4.0 * yCurrent) * (
    kr * (currentU + topU) - kl * (bottomU + currentU) +
    std::fabs(kr) * (currentU - topU) -
    std::fabs(kl) * (bottomU - currentU));

  Scalar const result = (1.0 - gamma) * secondOrder + gamma * firstOrder;

  return result;
}

template <typename Scalar>
inline Scalar
du2dx(Scalar const& currentU,
      Scalar const& leftU,
      Scalar const& rightU,
      Scalar const& currentX,
      Scalar const& leftX,
      Scalar const& rightX,
      Scalar const& gamma) {
  Scalar const dxShort = 0.5 * currentX;
  Scalar const dxLong0 = 0.5 * (leftX + currentX);
  Scalar const dxLong1 = 0.5 * (currentX + rightX);

  Scalar const kr = (dxLong1 - dxShort) / dxLong1 * currentU +
                    dxShort / dxLong1 * rightU;
  Scalar const kl = (dxLong0 - dxShort) / dxLong0 * currentU +
                    dxShort / dxLong0 * leftU;

  Scalar const secondOrder = (
    ((dxLong1 - dxShort) / dxLong1 * currentU +
     dxShort / dxLong1 * rightU) *
    ((dxLong1 - dxShort) / dxLong1 * currentU + dxShort / dxLong1 * rightU) -
    ((dxLong0 - dxShort) / dxLong0 * currentU + dxShort / dxLong0 * leftU) *
    ((dxLong0 - dxShort) / dxLong0 * currentU + dxShort / dxLong0 * leftU)
    ) / (2.0 * dxShort);

  Scalar const firstOrder = 1.0 / (4.0 * dxShort) *
                            (kr * (currentU + rightU) -
                             kl * (leftU + currentU) +
                             std::fabs(kr) * (currentU - rightU) -
                             std::fabs(kl) * (leftU - currentU));

  Scalar const result = (1.0 - gamma) * secondOrder + gamma * firstOrder;

  return result;
}

template <typename Scalar>
inline Scalar
computeFGH2D(
  Scalar const& currentU,
  Scalar const& currentV,
  Scalar const& leftU,
  Scalar const& rightU,
  Scalar const& rightV,
  Scalar const& bottomU,
  Scalar const& bottomV,
  Scalar const& topU,
  Scalar const& rightBottomV,
  Scalar const& currentX,
  Scalar const& currentY,
  Scalar const& leftX,
  Scalar const& rightX,
  Scalar const& bottomY,
  Scalar const& topY,
  Scalar const& re,
  Scalar const& gamma,
  Scalar const& gx,
  Scalar const& dt) {
  return currentU
         +  dt * (1.0 / re *
                  (d2udx2(currentU,
                          leftU,
                          rightU,
                          currentX,
                          rightX) +
                   d2udy2(currentU,
                          bottomU,
                          topU,
                          currentY,
                          bottomY,
                          topY)) -
                  du2dx(currentU,
                        leftU,
                        rightU,
                        currentX,
                        leftX,
                        rightX,
                        gamma) -
                  duvdy(currentU, // U(i, j, k)
                        currentV, // V(i, j, k)
                        bottomU, // U(i, j-1, k)
                        bottomV, // V(i, j-1, k)
                        rightBottomV, // V(i+1, j-1, k)
                        rightV, // V(i+1, j, k)
                        topU, // U(i, j+1, k)
                        currentY, // dy(i, j, k)
                        currentX, // dx(i, j, k)
                        bottomY, // dy(i, j-1, k)
                        topY, // dy(i, j+1, k)
                        rightX, // dx(i+1, j, k)
                        gamma) +
                  gx);
}

template <typename Scalar>
inline Scalar
computeFGH3D(
  Scalar const& currentU,
  Scalar const& currentV,
  Scalar const& currentW,
  Scalar const& leftU,
  Scalar const& rightU,
  Scalar const& rightV,
  Scalar const& rightW,
  Scalar const& bottomU,
  Scalar const& bottomV,
  Scalar const& topU,
  Scalar const& backU,
  Scalar const& backW,
  Scalar const& frontU,
  Scalar const& rightBottomV,
  Scalar const& rightBackW,
  Scalar const& currentX,
  Scalar const& currentY,
  Scalar const& currentZ,
  Scalar const& leftX,
  Scalar const& rightX,
  Scalar const& bottomY,
  Scalar const& topY,
  Scalar const& backZ,
  Scalar const& frontZ,
  Scalar const& re,
  Scalar const& gamma,
  Scalar const& gx,
  Scalar const& dt) {
  return currentU
         +  dt * (1.0 / re *
                  (d2udx2(currentU,
                          leftU,
                          rightU,
                          currentX,
                          rightX) +
                   d2udy2(currentU,
                          bottomU,
                          topU,
                          currentY,
                          bottomY,
                          topY) +
                   d2udy2(currentU, // d2udz2
                          backU,
                          frontU,
                          currentZ,
                          backZ,
                          frontZ)) -
                  du2dx(currentU,
                        leftU,
                        rightU,
                        currentX,
                        leftX,
                        rightX,
                        gamma) -
                  duvdy(currentU, // U(i, j, k)
                        currentV, // V(i, j, k)
                        bottomU, // U(i, j-1, k)
                        bottomV, // V(i, j-1, k)
                        rightBottomV, // V(i+1, j-1, k)
                        rightV, // V(i+1, j, k)
                        topU, // U(i, j+1, k)
                        currentY, // dy(i, j, k)
                        currentX, // dx(i, j, k)
                        bottomY, // dy(i, j-1, k)
                        topY, // dy(i, j+1, k)
                        rightX, // dx(i+1, j, k)
                        gamma) -
                  duvdy(currentU, // duwdz
                        currentW,
                        backU,
                        backW,
                        rightBackW,
                        rightW,
                        frontU,
                        currentZ,
                        currentX,
                        backZ,
                        frontZ,
                        rightX,
                        gamma) +
                  gx);
}

template <typename TCellAccessor,
          typename TSimulationParameters,
          typename Scalar,
          int D>
class FghProcessing {
  static inline void
  compute(TCellAccessor const&         accessor,
          TSimulationParameters const& simulationParameters,
          Scalar const&                dt) {}
};

template <typename TCellAccessor,
          typename TSimulationParameters,
          typename Scalar>
class FghProcessing<TCellAccessor, TSimulationParameters, Scalar, 3> {
  static inline void
  compute(TCellAccessor const&         accessor,
          TSimulationParameters const& simulationParameters,
          Scalar const&                dt) {
    auto current     = accessor.currentCell();
    auto left        = accessor.leftCell();
    auto right       = accessor.rightCell();
    auto bottom      = accessor.bottomCell();
    auto top         = accessor.topCell();
    auto back        = accessor.backCell();
    auto front       = accessor.currentCell();
    auto rightBottom = accessor.rightBottomCell();
    auto rightBack   = accessor.rightBackCell();
    auto topLeft     = accessor.topLeftCell();
    auto topBack     = accessor.topBackCell();
    auto leftFront   = accessor.leftFrontCell();
    auto bottomFront = accessor.bottomFrontCell();

    auto currentSize = accessor.currentSize();
    auto leftSize    = accessor.leftSize();
    auto rightSize   = accessor.rightSize();
    auto bottomSize  = accessor.bottomSize();
    auto topSize     = accessor.topSize();
    auto backSize    = accessor.backSize();
    auto frontSize   = accessor.current();

    current.fgh(0) =
      computeFGH3D(
        current->velocity(0),
        current->velocity(1),
        current->velocity(2),
        left->velocity(0),
        right->velocity(0),
        right->velocity(1),
        right->velocity(2),
        bottom->velocity(0),
        bottom->velocity(1),
        top->velocity(0),
        back->velocity(0),
        back->velocity(2),
        front->velocity(0),
        rightBottom->velocity(1),
        rightBack->velocity(2),
        currentSize(0),
        currentSize(1),
        currentSize(2),
        leftSize(0),
        rightSize(0),
        bottomSize(1),
        topSize(1),
        backSize(2),
        frontSize(2),
        simulationParameters.re(),
        simulationParameters.gamma(),
        simulationParameters.g(0),
        dt);

    current.fgh(1) =
      computeFGH3D(
        current->velocity(1),
        current->velocity(0),
        current->velocity(2),
        bottom->velocity(1),
        top->velocity(1),
        top->velocity(0),
        top->velocity(2),
        left->velocity(1),
        left->velocity(0),
        right->velocity(1),
        back->velocity(1),
        back->velocity(2),
        front->velocity(1),
        topLeft->velocity(0),
        topBack->velocity(2),
        currentSize(1),
        currentSize(0),
        currentSize(2),
        bottomSize(1),
        topSize(1),
        leftSize(0),
        rightSize(0),
        backSize(2),
        frontSize(2),
        simulationParameters.re(),
        simulationParameters.gamma(),
        simulationParameters.g(1),
        dt);

    current.fgh(2) =
      computeFGH3D(
        current->velocity(2),
        current->velocity(1),
        current->velocity(0),
        back->velocity(2),
        front->velocity(2),
        front->velocity(1),
        front->velocity(0),
        bottom->velocity(2),
        bottom->velocity(1),
        top->velocity(2),
        left->velocity(2),
        left->velocity(0),
        right->velocity(2),
        bottomFront->velocity(1),
        leftFront->velocity(0),
        currentSize(2),
        currentSize(1),
        currentSize(0),
        backSize(2),
        frontSize(2),
        bottomSize(1),
        topSize(1),
        leftSize(0),
        rightSize(0),
        simulationParameters.re(),
        simulationParameters.gamma(),
        simulationParameters.g(2),
        dt);
  }
};

template <typename TCellAccessor,
          typename TSimulationParameters,
          typename Scalar>
class FghProcessing<TCellAccessor, TSimulationParameters, Scalar, 2> {
  static inline void
  compute(TCellAccessor const&         accessor,
          TSimulationParameters const& simulationParameters,
          Scalar const&                dt) {
    auto current     = accessor.currentCell();
    auto left        = accessor.leftCell();
    auto right       = accessor.rightCell();
    auto bottom      = accessor.bottomCell();
    auto top         = accessor.topCell();
    auto rightBottom = accessor.rightBottomCell();
    auto topLeft     = accessor.topLeftCell();
    auto leftFront   = accessor.leftFrontCell();
    auto bottomFront = accessor.bottomFrontCell();

    auto currentSize = accessor.currentSize();
    auto leftSize    = accessor.leftSize();
    auto rightSize   = accessor.rightSize();
    auto bottomSize  = accessor.bottomSize();
    auto topSize     = accessor.topSize();

    current.fgh(0) =
      computeFGH2D(
        current->velocity(0),
        current->velocity(1),
        left->velocity(0),
        right->velocity(0),
        right->velocity(1),
        bottom->velocity(0),
        bottom->velocity(1),
        top->velocity(0),
        rightBottom->velocity(1),
        currentSize(0),
        currentSize(1),
        leftSize(0),
        rightSize(0),
        bottomSize(1),
        topSize(1),
        simulationParameters.re(),
        simulationParameters.gamma(),
        simulationParameters.g(0),
        dt);

    current.fgh(1) =
      computeFGH2D(
        current->velocity(1),
        current->velocity(0),
        bottom->velocity(1),
        top->velocity(1),
        top->velocity(0),
        left->velocity(1),
        left->velocity(0),
        right->velocity(1),
        topLeft->velocity(0),
        currentSize(1),
        currentSize(0),
        bottomSize(1),
        topSize(1),
        leftSize(0),
        rightSize(0),
        simulationParameters.re(),
        simulationParameters.gamma(),
        simulationParameters.g(1),
        dt);
  }
};
}

#endif
