#ifndef FsiSimulation_FluidSimulation_Configuration_hpp
#define FsiSimulation_FluidSimulation_Configuration_hpp

#include <Eigen/Core>

#include <array>
#include <memory>
#include <string>

namespace FsiSimulation {
namespace FluidSimulation {
class Configuration {
public:
  class Wall {
public:
    typedef double                      Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1> VectorDs;

public:
    Wall(int const& type)
      : _type(type) {}

    Wall(Wall const&) = delete;

    virtual
    ~Wall() {}

    Wall&
    operator=(Wall const&) = delete;

    bool
    operator==(int const& type) const {
      return _type == type;
    }

    int const&
    type() const {
      return _type;
    }

    virtual
    VectorDs const&
    velocity() const {
      static VectorDs const _velocity;

      return _velocity;
    }

    static int const Input          = 0;
    static int const ParabolicInput = 1;
    static int const Output         = 2;

private:
    int _type;
  };

  class Input : public Wall {
public:
    typedef Wall                    Base;
    typedef typename Base::Scalar   Scalar;
    typedef typename Base::VectorDs VectorDs;

public:
    Input(VectorDs const& velocity) : Base(Base::Input),
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
    ParabolicInput(VectorDs const& velocity) : Base(Base::ParabolicInput),
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
    Output() : Base(Base::Output) {}

    Output(Output const&) = delete;

    virtual
    ~Output() {}

    Output&
    operator=(Output const&) = delete;
  };

public:
  typedef double                          ScalarType;
  typedef Eigen::Matrix<int, 3, 1>        VectorDiType;
  typedef Eigen::Matrix<ScalarType, 3, 1> VectorDsType;
  typedef std::unique_ptr<Wall>           UniqueWallType;
  typedef std::array<std::array<UniqueWallType, 2>, 3>
    WallsType;

public:
  Configuration()
    : re(0),
      timeLimit(0),
      plotInterval(0),
      tau(0),
      gamma(0),
      dim(0) {}

  ScalarType   re;
  ScalarType   timeLimit;
  ScalarType   plotInterval;
  ScalarType   tau;
  ScalarType   gamma;
  int          dim;
  VectorDsType width;
  VectorDiType size;
  VectorDiType parallelizationSize;
  VectorDsType environment;
  WallsType    walls;
  std::string  filename;
};
}
}

#endif
