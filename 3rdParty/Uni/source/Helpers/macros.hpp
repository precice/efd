#pragma once

#define Uni_UNUSED(variable) \
  (void)(variable)

#define Uni_Property(type, name)                               \
  type const & name() const { return _ ## name; }              \
  type& name() { return _ ## name; }                           \
  void name(type const & name ## _) { _ ## name = name ## _; } \
private:                                                       \
  type _ ## name;

#define Uni_PublicProperty(type, name) \
  Uni_Property(type, name)             \
public:

#define Uni_PrivateProperty(type, name) \
public:                                 \
  Uni_Property(type, name)
