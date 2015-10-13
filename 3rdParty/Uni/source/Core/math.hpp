#ifndef Uni_Core_math_hpp
#define Uni_Core_math_hpp

namespace Uni {
template <typename T>
int
sign(T x) {
  return (T(0) < x) - (x < T(0));
}
}

#endif
