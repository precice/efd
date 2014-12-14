#ifndef FsiSimulation_FluidSimulation_Configuration_hpp
#define FsiSimulation_FluidSimulation_Configuration_hpp

#include <Eigen/Core>

#include <array>
#include <petscsys.h>
#include <string>

namespace FsiSimulation {
namespace FluidSimulation {
class TimestepParameters {
public:
  double dt;
  double tau;
};

class SmParameters {
public:
  double      finalTime;
  std::string type;
  std::string scenario;
};

class EnvironmentalParameters {
public:
  double gx;
  double gy;
  double gz;
};

class FlowParameters {
public:
  double Re;
};

class SolverParameters {
public:
  double gamma;
  int    maxIterations;
};

class GeometricParameters {
public:
  int dim;

  int sizeX;
  int sizeY;
  int sizeZ;

  double lengthX;
  double lengthY;
  double lengthZ;

  int meshsizeType;

  int stretchX;
  int stretchY;
  int stretchZ;
};

class WallParameters {
public:
  typedef Eigen::Matrix<double, 3, 1>            VectorDs;
  typedef std::array<std::array<VectorDs, 2>, 3> Velocities;
  typedef std::array<std::array<double, 2>, 3>   Pressures;

  double scalarLeft;
  double scalarRight;
  double scalarBottom;
  double scalarTop;
  double scalarFront;
  double scalarBack;

  Pressures pressures;

  double vectorLeft[3];
  double vectorRight[3];
  double vectorBottom[3];
  double vectorTop[3];
  double vectorFront[3];
  double vectorBack[3];

  Velocities velocities;
};

class VTKParameters {
public:
  double      interval;
  std::string prefix;
};

class StdOutParameters {
public:
  double interval;
};

class ParallelParameters {
public:
  int rank;

  int numProcessors[3];

  int leftNb;
  int rightNb;
  int bottomNb;
  int topNb;
  int frontNb;
  int backNb;

  int indices[3];
  int localSize[3];
  int firstCorner[3];

  PetscInt* sizes[3];
};

class BFStepParameters {
public:
  double xRatio;
  double yRatio;
};

class BaldwinLomaxModel {
public:
  double Uymax_up;
  double Umax_up;
  double ymax_up;
  double lmixmax_up;
  double Umax_down;
  double Uymax_down;
  double ymax_down;
  double lmixmax_down;
  double vel_tau;

  double      kappa;
  double      delta99;
  std::string modelType;
};

class Configuration {
public:
  SmParameters            simulation;
  TimestepParameters      timestep;
  EnvironmentalParameters environment;
  FlowParameters          flow;
  SolverParameters        solver;
  GeometricParameters     geometry;
  WallParameters          walls;
  VTKParameters           vtk;
  ParallelParameters      parallel;
  StdOutParameters        stdOut;
  BFStepParameters        bfStep;
  BaldwinLomaxModel       blm;
};
}
}

#endif
