#pragma once

#include <cstddef>

namespace FsiSimulation {
namespace FluidSimulation {
namespace Private {
template <typename TVecotrDi>
struct Hash {
  inline std::size_t
  operator()(TVecotrDi const& one) const {
    return one.sum();
  }
};
}
}
}
