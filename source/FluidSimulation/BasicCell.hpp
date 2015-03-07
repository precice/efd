#ifndef FsiSimulation_FluidSimulation_BasicCell_hpp
#define FsiSimulation_FluidSimulation_BasicCell_hpp

#include <Uni/Logging/macros>

#include <map>
#include <set>
#include <string>
#include <list>

namespace FsiSimulation {
namespace FluidSimulation {
class Attribute {
public:
  enum class Type {
    Vector = 0,
    Scalar = 1
  };

  Attribute(std::string name_,
            Type const& type_) :
    name(name_),
    type(type_) {}

  Attribute(std::string name_,
            int const&  index_,
            Type const& type_) :
    name(name_),
    index(index_),
    type(type_) {}

  Attribute(Attribute const& other) : name(other.name),
    index(other.index),
    type(other.type) {}

  std::string const name;
  int               index;
  Type const        type;
};

template <typename TDerived>
class BasicCellTraits {
public:
private:
  struct StringPointerLess {
    bool
    operator()(std::string const* string1, std::string const* string2) const {
      return *string1 < *string2;
    }
  };

  struct IntPointerLess {
    bool
    operator()(int const* int1, int const* int2) const {
      return *int1 < *int2;
    }
  };

  using AttributeIndexMap
          = std::map<int const*, Attribute*, IntPointerLess>;
  using AttributeNameMap
          = std::map<std::string const*, Attribute*, StringPointerLess>;
  using Attributes
          = std::list<Attribute>;

public:
  BasicCellTraits() {
    TDerived::initializeAttributes();
  }

  static void
  setAttribute(std::string const& name,
               int const&         index,
               Attribute::Type type) {
    _attributes.emplace_back(name, index, type);
    _attributeNameMap.insert(std::make_pair(
                                &_attributes.back().name,
                                &_attributes.back()));
    _attributeIndexMap.insert(std::make_pair(
                                 &_attributes.back().index,
                                 &_attributes.back()));
  }

  static int
  getAttributesSize() {
    return _attributes.size();
  }

  static Attribute*
  getAttribute(int const& index) {
    return _attributeIndexMap.find(&index)->second;
  }

  static Attribute*
  getAttribute(std::string const& name) {
    return _attributeNameMap.find(&name)->second;
  }

  static bool
  initializeAttributes() {
    static bool isInitialized = false;

    if (!isInitialized) {
      TDerived::initializeAttributes();
      isInitialized = true;
    }
  }

private:
  static Attributes        _attributes;
  static AttributeNameMap  _attributeNameMap;
  static AttributeIndexMap _attributeIndexMap;
};

template <typename TDerived>
typename BasicCellTraits<TDerived>::Attributes
BasicCellTraits<TDerived>::_attributes;

template <typename TDerived>
typename BasicCellTraits<TDerived>::AttributeNameMap
BasicCellTraits<TDerived>::_attributeNameMap;

template <typename TDerived>
typename BasicCellTraits<TDerived>::AttributeIndexMap
BasicCellTraits<TDerived>::_attributeIndexMap;

template <typename TCell>
class CellTraits {};
}
}

#endif
