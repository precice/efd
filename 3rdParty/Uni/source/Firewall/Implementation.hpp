#ifndef Uni_Firewall_Implementation_hpp
#define Uni_Firewall_Implementation_hpp

// Preprocessor {{{
// =============================================================================
#include "DeepConstPointer.hpp"

#define Uni_Firewall_IMPLEMENTATION_LINK_AS(Type, name) \
  friend class Type;                                    \
  Uni::Firewall::Implementation<Type> name;

#define Uni_Firewall_IMPLEMENTATION_LINK(Type) \
  Uni_Firewall_IMPLEMENTATION_LINK_AS(Type, _im)

#define Uni_Firewall_IMPLEMENTATION \
  class Private;                    \
  Uni_Firewall_IMPLEMENTATION_LINK(Private)
// =============================================================================
// }}} Preprocessor

namespace Uni {
namespace Firewall {
/**
 *  Provides a holder for an implementation pointer
 */
template <class T>
class Implementation : public DeepConstPointer<T> {
public:
  typedef T Interface;
public:
  Implementation() : DeepConstPointer<T>(new T) {}

  Implementation(T* pointer) : DeepConstPointer<T>(pointer) {}

  ~Implementation()
  { if (this->pointer()) { delete this->pointer(); } }

private:
  Implementation(Implementation const& other);

private:
  Implementation&
  operator=(Implementation const& other);
};
}
}

#endif
