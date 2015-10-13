#ifndef Uni_Firewall_DeepConstPointer_hpp
#define Uni_Firewall_DeepConstPointer_hpp

namespace Uni {
namespace Firewall {
/**
 *  Holds a pointer,
 *  Provides a correct return type of its pointer
 *  in accordance with constance of used scope.
 */
template <class T>
class DeepConstPointer {
public:
  typedef T Type;

public:
  DeepConstPointer(Type* pointer) : _pointer(pointer) {}

public:
  Type*
  pointer() { return _pointer; }

  Type const*
  pointer() const { return _pointer; }

  Type*
  operator->()
  { return _pointer; }

  Type const*
  operator->() const
  { return static_cast<Type const*>(_pointer); }

  Type&
  operator*()
  { return *_pointer; }

  Type const&
  operator*() const
  { return *_pointer; }

private:
  Type* _pointer;
};
}
}

#endif
