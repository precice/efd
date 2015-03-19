#ifndef FsiSimulation_FluidSimulation_Configuration_hpp
#define FsiSimulation_FluidSimulation_Configuration_hpp

#include <Eigen/Core>

#include <array>
#include <memory>
#include <string>

namespace FsiSimulation {
namespace FluidSimulation {
enum class SolverEnum {
  Sfsfd,
  Ifsfd,
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
    typedef double                      Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1> VectorDs;

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
  Configuration()
    : re(0),
    timeLimit(0),
    plotInterval(0),
    tau(0),
    gamma(0),
    dim(0),
    immersedBoundaryMethod(-1),
    solver(SolverEnum::Sfsfd) {
    for (int d = 0; d < 3; ++d) {
      walls[d][0] = UniqueWallType(new Input());
      walls[d][1] = UniqueWallType(new Input());
    }
  }

  ScalarType   re;
  ScalarType   timeLimit;
  ScalarType   plotInterval;
  int          iterationLimit;
  ScalarType   tau;
  ScalarType   gamma;
  int          dim;
  VectorDsType width;
  VectorDiType size;
  VectorDiType parallelizationSize;
  VectorDsType environment;
  WallTypes    walls;
  std::string  filename;
  int          immersedBoundaryMethod;
  float        alpha;
  SolverEnum   solver;

  static int const FeedbackForcingMethod;
};
}
}

#endif
