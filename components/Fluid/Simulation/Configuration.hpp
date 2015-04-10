#pragma once

#include <Uni/Logging/format>

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
  Configuration() {
    reset();
  }

  void
  reset() {
    re                  = 1;
    timeLimit           = 0;
    plotInterval        = 0;
    tau                 = 0.5;
    gamma               = 0.5;
    dimensions          = 2;
    width               = VectorDsType::Zero();
    size                = VectorDiType::Ones();
    parallelizationSize = VectorDiType::Ones();
    environment         = VectorDsType::Zero();
    solverType          = SolverEnum::Sfsfd;
    scalarType          = ScalarEnum::Double;
    outputType          = OutputEnum::Vtk;
    outerLayerSize      = 1;
    innerLayerSize      = 0;
    doImmersedBoundary  = false;
    doDebug             = false;

    for (int d = 0; d < 3; ++d) {
      walls[d][0].reset(new Input());
      walls[d][1].reset(new Input());
    }
  }

  std::string
  toString() {
    return Uni::Logging::format(
      "Re = {1}\n"
      "TimeLimit = {2}\n"
      "IterationLimit = {3}\n"
      "PlotInterval = {4}\n"
      "Tau = {5}\n"
      "Gamma = {6}\n"
      "Dimensions = {7}\n"
      "Width = {8}\n"
      "Size = {9}\n"
      "ParallelizationSize = {10}\n"
      "Environment = {11}\n"
      "Walls[0][0].velocity = {12}\n"
      "Walls[0][1].velocity = {13}\n"
      "Walls[1][0].velocity = {14}\n"
      "Walls[1][1].velocity = {15}\n"
      "Walls[2][0].velocity = {16}\n"
      "Walls[2][1].velocity = {17}\n",
      "OuterLayerSize = {18}\n",
      "InnerLayerSize = {19}\n",
      re,
      timeLimit,
      iterationLimit,
      plotInterval,
      tau,
      gamma,
      dimensions,
      width.transpose(),
      size.transpose(),
      parallelizationSize.transpose(),
      environment.transpose(),
      walls[0][0]->velocity().transpose(),
      walls[0][1]->velocity().transpose(),
      walls[1][0]->velocity().transpose(),
      walls[1][1]->velocity().transpose(),
      walls[2][0]->velocity().transpose(),
      walls[2][1]->velocity().transpose(),
      outerLayerSize,
      innerLayerSize);
  }

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
};
}
}
