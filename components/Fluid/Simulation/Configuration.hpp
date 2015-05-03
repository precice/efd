#pragma once

#include <Uni/ExecutionControl/exception>
#include <Uni/Firewall/Implementation>

#include <Eigen/Core>

#include <boost/any.hpp>

#include <array>
#include <memory>
#include <string>

namespace FsiSimulation {
namespace FluidSimulation {
class ConfigurationImplementation;

enum class SolverEnum {
  Sfsfd,
  Ifsfd,
};

enum class ScalarEnum {
  Float,
  Double,
  LongDouble
};

enum class OutputEnum {
  Vtk,
  Xdmf
};

enum class WallEnum {
  Input,
  ParabolicInput,
  Output
};
class Configuration {
public:
  class Wall {
public:
    using Scalar   = long double;
    using VectorDs = Eigen::Matrix<Scalar, 3, 1>;

public:
    Wall(WallEnum const& type)
      : _type(type) {}

    Wall(Wall const&) = delete;

    virtual
    ~Wall() {}

    Wall&
    operator=(Wall const&) = delete;

    bool
    operator==(WallEnum const& type) const {
      return _type == type;
    }

    WallEnum const&
    type() const {
      return _type;
    }

    virtual
    VectorDs const&
    velocity() const {
      static VectorDs const _velocity = VectorDs::Zero();

      return _velocity;
    }

private:
    WallEnum _type;
  };

  class Input : public Wall {
public:
    typedef Wall                    Base;
    typedef typename Base::Scalar   Scalar;
    typedef typename Base::VectorDs VectorDs;

public:
    Input() : Base(WallEnum::Input),
      _velocity(VectorDs::Zero()) {}

    Input(VectorDs const& velocity) : Base(WallEnum::Input),
      _velocity(velocity) {}

    Input(Input const&) = delete;

    virtual
    ~Input() {}

    Input&
    operator=(Input const&) = delete;

    VectorDs const&
    velocity() const {
      return _velocity;
    }

private:
    VectorDs _velocity;
  };

  class ParabolicInput : public Wall {
public:
    typedef Wall                    Base;
    typedef typename Base::Scalar   Scalar;
    typedef typename Base::VectorDs VectorDs;

public:
    ParabolicInput(VectorDs const& velocity) : Base(WallEnum::ParabolicInput),
      _velocity(velocity) {}

    ParabolicInput(ParabolicInput const&) = delete;

    virtual
    ~ParabolicInput() {}

    ParabolicInput&
    operator=(ParabolicInput const&) = delete;

    VectorDs const&
    velocity() const {
      return _velocity;
    }

private:
    VectorDs _velocity;
  };

  class Output : public Wall {
public:
    typedef Wall                    Base;
    typedef typename Base::Scalar   Scalar;
    typedef typename Base::VectorDs VectorDs;

public:
    Output() : Base(WallEnum::Output) {}

    Output(Output const&) = delete;

    virtual
    ~Output() {}

    Output&
    operator=(Output const&) = delete;
  };

public:
  typedef double                                       ScalarType;
  typedef Eigen::Matrix<int, 3, 1>                     VectorDiType;
  typedef Eigen::Matrix<ScalarType, 3, 1>              VectorDsType;
  typedef std::unique_ptr<Wall>                        UniqueWallType;
  typedef std::array<std::array<UniqueWallType, 2>, 3> WallTypes;

public:
  Configuration();

  ~Configuration();

  void
  reset();

  std::string
  toString() const;

  ScalarType   re;
  ScalarType   timeLimit;
  ScalarType   plotInterval;
  int          iterationLimit;
  ScalarType   tau;
  ScalarType   gamma;
  int          dimensions;
  VectorDsType width;
  VectorDiType size;
  VectorDiType parallelizationSize;
  VectorDsType environment;
  WallTypes    walls;
  std::string  filename;
  SolverEnum   solverType;
  ScalarEnum   scalarType;
  OutputEnum   outputType;
  unsigned     outerLayerSize;
  unsigned     innerLayerSize;
  bool         doImmersedBoundary;
  bool         doDebug;

  bool
  is(std::string const& name) const;

  template <typename T>
  T
  get(std::string const& name) const {
    try {
      return boost::any_cast<T>(_find(name));
    } catch (boost::bad_any_cast const&) {
      logFatal("Cannot retrieve property '{1}' value", name);
      throw std::exception();
    }
  }

  template <typename T>
  T
  tryToGet(std::string const& name) const {
    return boost::any_cast<T>(_find(name));
  }

  template <typename T>
  bool
  isOfType(std::string const& name) const {
    auto value = _find(name);

    try {
      boost::any_cast<T>(value);

      return true;
    } catch (boost::bad_any_cast const&) {
      return false;
    }
  }

  void
  set(std::string const& name, boost::any const& value);

  void
  set(std::string const& name);

private:
  void
  _set(std::string const& name, boost::any const& value);

  boost::any const&
  _find(std::string const& name) const;

  boost::any&
  _find(std::string const& name);

  Uni_Firewall_IMPLEMENTATION_LINK(ConfigurationImplementation);
};
}
}
