#ifndef Fluid_Simulation_function_hpp
#define Fluid_Simulation_function_hpp

#include <Eigen/Core>

#include <Uni/Logging/macros>

#include <cmath>
#include <limits>

namespace Fluid {
namespace Simulation {

template <typename TScalar>
inline TScalar
dudx(TScalar const& leftU,
     TScalar const& rightU,
     TScalar const& currentX,
     TScalar const& leftX,
     TScalar const& rightX) {
  TScalar const width = (leftX + rightX) / 2.0;

  return (rightU - leftU) / (currentX + width);
}

template <typename TScalar>
inline TScalar
d2udx2(TScalar const& currentU,
       TScalar const& leftU,
       TScalar const& rightU,
       TScalar const& currentX,
       TScalar const& rightX) {
  TScalar const factor = currentX + rightX;

  return 2.0 * (rightU / (rightX * factor) -
                currentU / (currentX * rightX) +
                leftU / (currentX * factor));
}

template <typename TScalar>
inline TScalar
d2udy2(TScalar const& currentU,
       TScalar const& bottomU,
       TScalar const& topU,
       TScalar const& currentY,
       TScalar const& bottomY,
       TScalar const& topY) {
  TScalar const factorBottom = 0.5 * (currentY + bottomY);
  TScalar const factorTop    = 0.5 * (currentY + topY);
  TScalar const factor       = factorBottom + factorTop;

  return 2.0 * (topU / (factorTop * factor) -
                currentU / (factorBottom * factorTop) +
                bottomU / (factorBottom * factor));
}

template <typename TScalar>
inline TScalar
duvdy_2(TScalar const& currentU,
        TScalar const& currentV,
        TScalar const& bottomU,
        TScalar const& bottomV,
        TScalar const& rightBottomV,
        TScalar const& rightV,
        TScalar const& topU,
        TScalar const& currentY,
        TScalar const& currentX,
        TScalar const& bottomY,
        TScalar const& topY,
        TScalar const& rightX,
        TScalar const& gamma) {
  TScalar const xCurrent = 0.5 * currentX;
  TScalar const yCurrent = 0.5 * currentY;
  TScalar const yBottom  = 0.5 * (currentY + bottomY);
  TScalar const yTop     = 0.5 * (currentY + topY);
  TScalar const xRight   = 0.5 * (currentX + rightX);

  TScalar const ar = (xRight - xCurrent) / xRight * currentV
                     + xCurrent / xRight * rightV;
  TScalar const al = (xRight - xCurrent) / xRight * bottomV
                     + xCurrent / xRight * rightBottomV;

  TScalar const br = (yTop - yCurrent) / yTop * currentU
                     + yCurrent / yTop * topU;
  TScalar const bl = (yBottom - yCurrent) / yBottom * currentU
                     + yCurrent / yBottom * bottomU;

  TScalar const cr = (yTop - yCurrent) / yTop * currentU
                     - yCurrent / yTop * topU;
  TScalar const cl = yCurrent / yBottom * bottomU
                     - (yBottom - yCurrent) / yBottom * currentU;

  TScalar const result
    = ((ar * br - al * bl)
       + gamma * (std::abs(ar) * cr - std::abs(al) * cl))
      / (2.0 * yCurrent);

  return result;
}

template <typename TScalar>
inline TScalar
duvdy(TScalar const& currentU,
      TScalar const& currentV,
      TScalar const& bottomU,
      TScalar const& bottomV,
      TScalar const& rightBottomV,
      TScalar const& rightV,
      TScalar const& topU,
      TScalar const& currentY,
      TScalar const& currentX,
      TScalar const& bottomY,
      TScalar const& topY,
      TScalar const& rightX,
      TScalar const& gamma) {
  TScalar const yCurrent = 0.5 * currentY;
  TScalar const yBottom  = 0.5 * (currentY + bottomY);
  TScalar const yTop     = 0.5 * (currentY + topY);
  TScalar const xCurrent = 0.5 * currentX;
  TScalar const xRight   = 0.5 * (currentX + rightX);

  TScalar const kr = (xRight - xCurrent) / xRight * currentV
                     + xCurrent / xRight * rightV;
  TScalar const kl = (xRight - xCurrent) / xRight * bottomV
                     + xCurrent / xRight * rightBottomV;

  TScalar const secondOrder
    = (kr * ((yTop - yCurrent) / yTop * currentU
             + yCurrent / yTop * topU)
       - kl * ((yBottom - yCurrent) / yBottom * currentU
               + yCurrent / yBottom * bottomU)
       ) / (2.0 * yCurrent);

  TScalar const firstOrder
    = (kr * (currentU + topU) - kl * (bottomU + currentU)
       + std::fabs(kr) * (currentU - topU)
       - std::fabs(kl) * (bottomU - currentU)) / (4.0 * yCurrent);

  TScalar const result = (1.0 - gamma) * secondOrder + gamma * firstOrder;

  return result;
}

template <typename TScalar>
inline TScalar
du2dx(TScalar const& currentU,
      TScalar const& leftU,
      TScalar const& rightU,
      TScalar const& currentX,
      TScalar const& leftX,
      TScalar const& rightX,
      TScalar const& gamma) {
  TScalar const dxShort = 0.5 * currentX;
  TScalar const dxLong0 = 0.5 * (leftX + currentX);
  TScalar const dxLong1 = 0.5 * (currentX + rightX);

  TScalar const kr
    = (dxLong1 - dxShort) / dxLong1 * currentU + dxShort / dxLong1 * rightU;

  TScalar const kl
    = (dxLong0 - dxShort) / dxLong0 * currentU + dxShort / dxLong0 * leftU;

  TScalar const secondOrder = (kr * kr - kl * kl) / (2.0 * dxShort);

  TScalar const firstOrder
    = (kr * (currentU + rightU) -
       kl * (leftU + currentU) +
       std::fabs(kr) * (currentU - rightU) -
       std::fabs(kl) * (leftU - currentU)) / (4.0 * dxShort);

  TScalar const result = (1.0 - gamma) * secondOrder + gamma * firstOrder;

  return result;
}

template <typename TScalar>
inline TScalar
computeDiffusion2D(TScalar const& currentU,
                   TScalar const& leftU,
                   TScalar const& rightU,
                   TScalar const& bottomU,
                   TScalar const& topU,
                   TScalar const& currentX,
                   TScalar const& currentY,
                   TScalar const& rightX,
                   TScalar const& bottomY,
                   TScalar const& topY) {
  return d2udx2(currentU,
                leftU,
                rightU,
                currentX,
                rightX)
         + d2udy2(currentU,
                  bottomU,
                  topU,
                  currentY,
                  bottomY,
                  topY);
}

template <typename TScalar>
inline TScalar
computeDiffusion3D(TScalar const& currentU,
                   TScalar const& leftU,
                   TScalar const& rightU,
                   TScalar const& bottomU,
                   TScalar const& topU,
                   TScalar const& backU,
                   TScalar const& frontU,
                   TScalar const& currentX,
                   TScalar const& currentY,
                   TScalar const& currentZ,
                   TScalar const& rightX,
                   TScalar const& bottomY,
                   TScalar const& topY,
                   TScalar const& backZ,
                   TScalar const& frontZ) {
  return d2udx2(currentU,
                leftU,
                rightU,
                currentX,
                rightX)
         + d2udy2(currentU,
                  bottomU,
                  topU,
                  currentY,
                  bottomY,
                  topY)
         + d2udy2(currentU, // d2udz2
                  backU,
                  frontU,
                  currentZ,
                  backZ,
                  frontZ);
}

template <typename TScalar>
inline TScalar
computeConvection2D(TScalar const& currentU,
                    TScalar const& currentV,
                    TScalar const& leftU,
                    TScalar const& rightU,
                    TScalar const& rightV,
                    TScalar const& bottomU,
                    TScalar const& bottomV,
                    TScalar const& topU,
                    TScalar const& rightBottomV,
                    TScalar const& currentX,
                    TScalar const& currentY,
                    TScalar const& leftX,
                    TScalar const& rightX,
                    TScalar const& bottomY,
                    TScalar const& topY,
                    TScalar const& gamma) {
  return du2dx(currentU,
               leftU,
               rightU,
               currentX,
               leftX,
               rightX,
               gamma)
         + duvdy(currentU,
                 currentV,
                 bottomU,
                 bottomV,
                 rightBottomV,
                 rightV,
                 topU,
                 currentY,
                 currentX,
                 bottomY,
                 topY,
                 rightX,
                 gamma);
}

template <typename TScalar>
inline TScalar
computeConvection3D(TScalar const& currentU,
                    TScalar const& currentV,
                    TScalar const& currentW,
                    TScalar const& leftU,
                    TScalar const& rightU,
                    TScalar const& rightV,
                    TScalar const& rightW,
                    TScalar const& bottomU,
                    TScalar const& bottomV,
                    TScalar const& topU,
                    TScalar const& backU,
                    TScalar const& backW,
                    TScalar const& frontU,
                    TScalar const& rightBottomV,
                    TScalar const& rightBackW,
                    TScalar const& currentX,
                    TScalar const& currentY,
                    TScalar const& currentZ,
                    TScalar const& leftX,
                    TScalar const& rightX,
                    TScalar const& bottomY,
                    TScalar const& topY,
                    TScalar const& backZ,
                    TScalar const& frontZ,
                    TScalar const& gamma) {
  return du2dx(currentU,
               leftU,
               rightU,
               currentX,
               leftX,
               rightX,
               gamma)
         + duvdy(currentU,
                 currentV,
                 bottomU,
                 bottomV,
                 rightBottomV,
                 rightV,
                 topU,
                 currentY,
                 currentX,
                 bottomY,
                 topY,
                 rightX,
                 gamma)
         + duvdy(currentU, // duwdz
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
                 gamma);
}

template <int TD>
class DiffusionProcessing {
  template <typename TCellAccessor>
  static inline typename TCellAccessor::VectorDsType
  compute(TCellAccessor const&) {}
};

template <>
class DiffusionProcessing<2> {
public:
  template <typename TCellAccessor>
  static inline
  typename TCellAccessor::VectorDsType
  compute(TCellAccessor const& accessor) {
    typedef typename TCellAccessor::VectorDsType Vector;
    Vector result;

    for (int d1 = 0; d1 < 2; ++d1) {
      int d2 = d1 + 1;

      if (d1 == 1) {
        d2 = 0;
      }

      result(d1) = computeDiffusion2D(
        accessor.velocity(d1),
        accessor.velocity(d1, -1, d1),
        accessor.velocity(d1, +1, d1),
        accessor.velocity(d2, -1, d1),
        accessor.velocity(d2, +1, d1),
        accessor.width(d1),
        accessor.width(d2),
        accessor.width(d1, +1, d1),
        accessor.width(d2, -1, d2),
        accessor.width(d2, +1, d2));
    }

    return result;
  }
};

template <>
class DiffusionProcessing<3> {
public:
  template <typename TCellAccessor>
  static inline
  typename TCellAccessor::VectorDsType
  compute(TCellAccessor const& accessor) {
    typedef typename TCellAccessor::VectorDsType Vector;

    Vector result;

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

      result(d1) = computeDiffusion3D(
        accessor.velocity(d1),
        accessor.velocity(d1, -1, d1),
        accessor.velocity(d1, +1, d1),
        accessor.velocity(d2, -1, d1),
        accessor.velocity(d2, +1, d1),
        accessor.velocity(d3, -1, d1),
        accessor.velocity(d3, +1, d1),
        accessor.width(d1),
        accessor.width(d2),
        accessor.width(d3),
        accessor.width(d1, +1, d1),
        accessor.width(d2, -1, d2),
        accessor.width(d2, +1, d2),
        accessor.width(d3, -1, d3),
        accessor.width(d3, +1, d3));
    }

    return result;
  }
};

template <int TD>
class ConvectionProcessing {
  template <typename TCellAccessor>
  static inline typename TCellAccessor::VectorDsType
  compute(TCellAccessor const&,
          typename TCellAccessor::ScalarType const&) {}
};

template <>
class ConvectionProcessing<2> {
public:
  template <typename TCellAccessor>
  static inline typename TCellAccessor::VectorDsType
  compute(TCellAccessor const&                      accessor,
          typename TCellAccessor::ScalarType const& gamma) {
    typedef typename TCellAccessor::VectorDsType Vector;

    Vector result;

    for (int d1 = 0; d1 < 2; ++d1) {
      int d2 = d1 + 1;

      if (d1 == 1) {
        d2 = 0;
      }
      result(d1) = computeConvection2D(
        accessor.velocity(d1),
        accessor.velocity(d2),
        accessor.velocity(d1, -1, d1),
        accessor.velocity(d1, +1, d1),
        accessor.velocity(d1, +1, d2),
        accessor.velocity(d2, -1, d1),
        accessor.velocity(d2, -1, d2),
        accessor.velocity(d2, +1, d1),
        accessor.velocity(d2, -1, d1, +1, d2),
        accessor.width(d1),
        accessor.width(d2),
        accessor.width(d1, -1, d1),
        accessor.width(d1, +1, d1),
        accessor.width(d2, -1, d2),
        accessor.width(d2, +1, d2),
        gamma);
    }

    return result;
  }
};

template <>
class ConvectionProcessing<3> {
public:
  template <typename TCellAccessor>
  static inline typename TCellAccessor::VectorDsType
  compute(TCellAccessor const&                      accessor,
          typename TCellAccessor::ScalarType const& gamma) {
    typedef typename TCellAccessor::VectorDsType Vector;
    Vector result;

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

      result(d1) = computeConvection3D(
        accessor.velocity(d1),
        accessor.velocity(d2),
        accessor.velocity(d3),
        accessor.velocity(d1, -1, d1),
        accessor.velocity(d1, +1, d1),
        accessor.velocity(d1, +1, d2),
        accessor.velocity(d1, +1, d3),
        accessor.velocity(d2, -1, d1),
        accessor.velocity(d2, -1, d2),
        accessor.velocity(d2, +1, d1),
        accessor.velocity(d3, -1, d1),
        accessor.velocity(d3, -1, d3),
        accessor.velocity(d3, +1, d1),
        accessor.velocity(d2, -1, d1, +1, d2),
        accessor.velocity(d3, -1, d1, +1, d3),
        accessor.width(d1),
        accessor.width(d2),
        accessor.width(d3),
        accessor.width(d1, -1, d1),
        accessor.width(d1, +1, d1),
        accessor.width(d2, -1, d2),
        accessor.width(d2, +1, d2),
        accessor.width(d3, -1, d3),
        accessor.width(d3, +1, d3),
        gamma);
    }

    return result;
  }
};

template <int TSolverType>
struct PressureProcessing {};

template <>
struct PressureProcessing<0> {
  template <typename TCellAccessor>
  static inline typename TCellAccessor::VectorDsType
  grad(TCellAccessor const& accessor) {
    typedef typename TCellAccessor::VectorDsType Vector;
    Vector result;

    for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
      result(d)
        = (accessor.pressure(d, +1) - accessor.pressure())
          / (0.5 * (accessor.width(d, +1, d) + accessor.width(d)));
    }

    return result;
  }

  template <typename TCellAccessor, typename TScalar>
  static void
  initialize(TCellAccessor const& accessor, TScalar const& value) {
    accessor.pressure() = value;
  }
};

template <>
struct PressureProcessing<1> {
  template <typename TCellAccessor>
  static inline typename TCellAccessor::VectorDsType
  grad(TCellAccessor const& accessor) {
    typedef typename TCellAccessor::VectorDsType Vector;
    Vector result;

    for (int d = 0; d < TCellAccessor::Dimensions; ++d) {
      result(d)
        = (accessor.projectionTerm(d, +1) - accessor.projectionTerm())
          / (0.5 * (accessor.width(d, +1, d) + accessor.width(d)));
    }

    return result;
  }

  template <typename TCellAccessor, typename TScalar>
  static void
  initialize(TCellAccessor const& accessor, TScalar const& value) {
    accessor.pressure()       = value;
    accessor.projectionTerm() = value;
  }
};
}
}
#endif
