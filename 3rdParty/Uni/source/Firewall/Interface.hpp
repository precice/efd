#ifndef Uni_Firewall_Interface_hpp
#define Uni_Firewall_Interface_hpp

// Preprocessor {{{
// =============================================================================
#include "DeepConstPointer.hpp"

#define Uni_Firewall_INTERFACE_DECLARE(Type) \
  friend class Type;

#define Uni_Firewall_INTERFACE_LINK_AS(Type, name) \
  friend class Type;                               \
  Uni::Firewall::Interface<Type> name;

#define Uni_Firewall_INTERFACE_LINK(Type) \
  Uni_Firewall_INTERFACE_LINK_AS(Type, _in)
// =============================================================================
// }}} Preprocessor

namespace Uni {
namespace Firewall {
template <class T>
class Interface : public DeepConstPointer<T> {
public:
  Interface(T* pointer) : DeepConstPointer<T>(pointer) {}

private:
  Interface(Interface const& other);

private:
  Interface&
  operator=(Interface const& other);
};
}
}

#endif
