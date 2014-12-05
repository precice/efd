#ifndef FsiSimulation_mystencils_hpp
#define FsiSimulation_mystencils_hpp

#include <Eigen/Core>

#include <Uni/Logging/macros>

#include <mpi.h>
#include <petsc.h>

#include <cmath>
#include <limits>

namespace FsiSimulation {
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

template <typename TCellAccessor, typename Scalar, int D>
void
computeMaxVelocity(TCellAccessor const&                    accessor,
                   typename TCellAccessor::Cell::Velocity& maxVelocity) {
  maxVelocity =
    accessor.currentCell()->velocity().cwiseAbs().cwiseQuotient(
      accessor.currentWidth()).cwiseMax(maxVelocity);

  // if (accessor.currentCell()->velocity().maxCoeff() > 27.0) {
  // logInfo("Max {1}\n{2}\n{3}",
  // accessor.indexValues().transpose(),
  // accessor.currentCell()->velocity().transpose(),
  // accessor.currentCell()->velocity().cwiseAbs().transpose());
  // }
}

template <typename Scalar, int D>
struct TimeStepProcessing {
  typedef  Eigen::Matrix<Scalar, D, 1> VectorDs;
  static inline Scalar
  compute(Scalar const&   re,
          Scalar const&   tau,
          VectorDs const& minCellWidth,
          VectorDs const& maxVelocity) {
    Scalar localMin, globalMin;
    Scalar factor = 0.0;

    factor = minCellWidth.cwiseProduct(minCellWidth).cwiseInverse().sum();

    localMin = std::min((Scalar)(re / (2.0 * factor)),
                        maxVelocity.cwiseInverse().minCoeff());

    globalMin = std::numeric_limits<Scalar>::max();
    mpiAllReduceMin(localMin, globalMin);

    factor  = globalMin;
    factor *= tau;

    return factor;
  }
};

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

  Scalar const secondOrder =
    (((xRight - xCurrent) / xRight * currentV + xCurrent / xRight * rightV) *
     ((yTop - yCurrent) / yTop * currentU + yCurrent / yTop * topU) -
     ((xRight - xCurrent) / xRight * bottomV + xCurrent / xRight *
      rightBottomV) *
     ((yBottom - yCurrent) / yBottom * currentU + yCurrent / yBottom * bottomU)
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

  Scalar const kr =
    (dxLong1 - dxShort) / dxLong1 * currentU + dxShort / dxLong1 * rightU;

  Scalar const kl =
    (dxLong0 - dxShort) / dxLong0 * currentU + dxShort / dxLong0 * leftU;

  Scalar const secondOrder = (kr * kr - kl * kl) / (2.0 * dxShort);

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
computeFGH2D(Scalar const& currentU,
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
                          rightX)
                   + d2udy2(currentU,
                            bottomU,
                            topU,
                            currentY,
                            bottomY,
                            topY))
                  - du2dx(currentU,
                          leftU,
                          rightU,
                          currentX,
                          leftX,
                          rightX,
                          gamma)
                  - duvdy(currentU, // U(i, j, k)
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
                          gamma)
                  + gx);
}

template <typename Scalar>
inline Scalar
computeFGH3D(Scalar const& currentU,
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
public:
  static inline void
  compute(TCellAccessor const&         accessor,
          TSimulationParameters const& simulationParameters,
          Scalar const&                dt) {
    for (int d1 = 0; d1 < 3; ++d1) {
      int d2 = d1 + 1;
      int d3 = d2 + 1;

      if (d1 == 1) {
        d2 = 0;
        d3 = 2;
      } else if (d1 == 2) {
        d2 = 0;
        d3 = 1;
      }
      accessor.currentCell()->fgh(d1) = computeFGH3D(
        accessor.currentCell()->velocity(d1),
        accessor.currentCell()->velocity(d2),
        accessor.currentCell()->velocity(d3),
        accessor.leftCellInDimension(d1)->velocity(d1),
        accessor.rightCellInDimension(d1)->velocity(d1),
        accessor.rightCellInDimension(d1)->velocity(d2),
        accessor.rightCellInDimension(d1)->velocity(d3),
        accessor.leftCellInDimension(d2)->velocity(d1),
        accessor.leftCellInDimension(d2)->velocity(d2),
        accessor.rightCellInDimension(d2)->velocity(d1),
        accessor.leftCellInDimension(d3)->velocity(d1),
        accessor.leftCellInDimension(d3)->velocity(d3),
        accessor.rightCellInDimension(d3)->velocity(d1),
        accessor.leftRightCellInDimensions(d2, d1)->velocity(d2),
        accessor.leftRightCellInDimensions(d3, d1)->velocity(d3),
        accessor.currentWidth() (d1),
        accessor.currentWidth() (d2),
        accessor.currentWidth() (d3),
        accessor.leftWidthInDimension (d1)(d1),
        accessor.rightWidthInDimension (d1)(d1),
        accessor.leftWidthInDimension (d2)(d2),
        accessor.rightWidthInDimension (d2)(d2),
        accessor.leftWidthInDimension (d3)(d3),
        accessor.rightWidthInDimension (d3)(d3),
        simulationParameters.re(),
        simulationParameters.gamma(),
        simulationParameters.g(d1),
        dt);
    }
  }
};

template <typename TCellAccessor,
          typename TSimulationParameters,
          typename Scalar>
class FghProcessing<TCellAccessor, TSimulationParameters, Scalar, 2> {
public:
  static inline void
  compute(TCellAccessor const&         accessor,
          TSimulationParameters const& simulationParameters,
          Scalar const&                dt) {
    for (int d1 = 0; d1 < 2; ++d1) {
      int d2 = d1 + 1;

      if (d1 == 1) {
        d2 = 0;
      }
      accessor.currentCell()->fgh(d1) = computeFGH2D(
        accessor.currentCell()->velocity(d1),
        accessor.currentCell()->velocity(d2),
        accessor.leftCellInDimension(d1)->velocity(d1),
        accessor.rightCellInDimension(d1)->velocity(d1),
        accessor.rightCellInDimension(d1)->velocity(d2),
        accessor.leftCellInDimension(d2)->velocity(d1),
        accessor.leftCellInDimension(d2)->velocity(d2),
        accessor.rightCellInDimension(d2)->velocity(d1),
        accessor.leftRightCellInDimensions(d2, d1)->velocity(d2),
        accessor.currentWidth() (d1),
        accessor.currentWidth() (d2),
        accessor.leftWidthInDimension (d1)(d1),
        accessor.rightWidthInDimension (d1)(d1),
        accessor.leftWidthInDimension (d2)(d2),
        accessor.rightWidthInDimension (d2)(d2),
        simulationParameters.re(),
        simulationParameters.gamma(),
        simulationParameters.g(d1),
        dt);
    }
  }
};

template <typename TCellAccessor, typename Scalar, int D>
class RhsProcessing {
public:
  static inline Scalar
  compute(TCellAccessor const& accessor,
          Scalar const&        dt) {
    auto temp(accessor.currentCell()->fgh());

    for (int d = 0; d < D; ++d) {
      temp(d) -= accessor.leftCellInDimension(d)->fgh(d);
    }
    temp = temp.cwiseQuotient(accessor.currentWidth());

    Scalar result = temp.sum();
    result = 1.0 / dt * result;

    return result;
  }
};

template <typename TGrid, typename TCellAccessor, typename Scalar, int D>
class PressurePoissonStencilProcessing {
public:
  template <typename TStencil, typename TRow, typename TColumns>
  static inline void
  compute(TGrid const&         grid,
          TCellAccessor const& accessor,
          TStencil&            stencil,
          TRow&                row,
          TColumns&            columns) {
    typedef Eigen::Matrix<Scalar, 2* D, 1> Vector2Ds;

    Vector2Ds meanWidths;

    for (int d = 0; d < D; ++d) {
      meanWidths(2 * d) = 0.5 *
                          (accessor.currentWidth() (d) +
                           accessor.leftWidthInDimension (d)(d));
      meanWidths(2 * d + 1) = 0.5 *
                              (accessor.currentWidth() (d) +
                               accessor.rightWidthInDimension (d)(d));

      auto leftIndex(accessor.leftIndexInDimension(d)); // - grid.leftIndent());
      auto rightIndex(accessor.rightIndexInDimension(d)); // -
                                                          // grid.leftIndent());
      columns[2 * d].i     = leftIndex(0);
      columns[2 * d].j     = leftIndex(1);
      columns[2 * d + 1].i = rightIndex(0);
      columns[2 * d + 1].j = rightIndex(1);

      if (D == 3) {
        columns[2 * d].k     = leftIndex(2);
        columns[2 * d + 1].k = rightIndex(2);
      }
    }
    auto currentIndex(accessor.currentIndex());
    // currentIndex    -= grid.leftIndent();
    columns[2 * D].i = currentIndex(0);
    columns[2 * D].j = currentIndex(1);

    if (D == 3) {
      columns[2 * D].k = currentIndex(2);
    }
    row = columns[2 * D];

    stencil[2 * D] = 0;

    for (int d = 0; d < D; ++d) {
      auto meanLeftAndRightWidth = meanWidths(2 * d) + meanWidths(2 * d + 1);
      stencil[2 * d] =
        2.0 / (meanWidths(2 * d) * meanLeftAndRightWidth);
      stencil[2 * d + 1] =
        2.0 / (meanWidths(2 * d + 1) * meanLeftAndRightWidth);
      stencil[2 * D] -= 2.0 / (meanWidths(2 * d) * meanWidths(2 * d + 1));
    }
  }
};

template <typename TCellAccessor, typename Scalar, int D>
class VelocityProcessing {
public:
  static inline void
  compute(TCellAccessor const& accessor,
          Scalar const&        dt) {
    for (int d = 0; d < D; ++d) {
      accessor.currentCell()->velocity() (d) =
        accessor.currentCell()->fgh() (d) -
        dt / (0.5 * (accessor.rightWidthInDimension (d)(d) +
                     accessor.currentWidth() (d)))
        * (accessor.rightCellInDimension(d)->pressure() -
           accessor.currentCell()->pressure());
    }
  }
};
}

#endif
