#include "Configuration.hpp"

#include "TinyXml2/tinyxml2.h"

#include <Uni/Logging/macros>

using FsiSimulation::EntryPoint::Configuration;

void
readFloatMandatory(double&               storage,
                   tinyxml2::XMLElement* node,
                   char const*           tag) {
  double value;

  if (node->QueryDoubleAttribute(tag, &value) != tinyxml2::XML_NO_ERROR) {
    logError("Error while reading mandatory argument");
  } else {
    storage = (double)value;
  }
}

void
readFloatOptional(double&               storage,
                  tinyxml2::XMLElement* node,
                  char const*           tag,
                  double                defaultValue = 0) {
  double value;
  int    result = node->QueryDoubleAttribute(tag, &value);

  if (result == tinyxml2::XML_NO_ATTRIBUTE) {
    storage = defaultValue;
  } else if (result == tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
    logError("Error while reading optional argument");
  } else {
    storage = (double)value;
  }
}

void
readIntMandatory(int&                  storage,
                 tinyxml2::XMLElement* node,
                 char const*           tag) {
  int value;

  if (node->QueryIntAttribute(tag, &value) != tinyxml2::XML_NO_ERROR) {
    logError("Error while reading mandatory argument");
  } else {
    storage = value;
  }
}

void
readIntOptional(int&                  storage,
                tinyxml2::XMLElement* node,
                char const*           tag,
                int                   defaultValue = 0) {
  int value;
  int result = node->QueryIntAttribute(tag, &value);

  storage = value;

  if (result == tinyxml2::XML_NO_ATTRIBUTE) {
    storage = defaultValue;
  } else if (result == tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
    logError("Error while reading optional argument");
  }
}

void
readBoolMandatory(bool&                 storage,
                  tinyxml2::XMLElement* node,
                  char const*           tag) {
  bool value;

  if (node->QueryBoolAttribute(tag, &value) != tinyxml2::XML_NO_ERROR) {
    logError("Error while reading mandatory argument");
  } else {
    storage = value;
  }
}

void
readBoolOptional(bool&                 storage,
                 tinyxml2::XMLElement* node,
                 char const*           tag,
                 bool                  defaultValue = false) {
  int result = node->QueryBoolAttribute(tag, &storage);

  if (result == tinyxml2::XML_NO_ATTRIBUTE) {
    storage = defaultValue;
  } else if (result == tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
    logError("Error while reading optional argument");
  }
}

void
readStringMandatory(std::string&          storage,
                    tinyxml2::XMLElement* node) {
  char const* myText = node->GetText();

  if (myText == NULL) {
    std::string const nodename = node->Name();
    std::cerr << "ERROR in file " << __FILE__ << ", line " << __LINE__ << ": ";
    std::cerr << "No string specified for this node: " << nodename << std::endl;
    exit(2);
  } else {
    storage = node->GetText();

    if (!storage.compare("")) {
      logError("Missing mandatory string!");
    }
  }
}

void
readWall(tinyxml2::XMLElement* wall,
         double*               vector,
         double&               scalar) {
  tinyxml2::XMLElement* quantity = wall->FirstChildElement("vector");

  if (quantity != NULL) {
    readFloatOptional(vector[0], quantity, "x");
    readFloatOptional(vector[1], quantity, "y");
    readFloatOptional(vector[2], quantity, "z");
  }
  quantity = wall->FirstChildElement("scalar");

  if (quantity != NULL) {
    readFloatOptional(scalar, quantity, "value");
  }
}

void
broadcastString(std::string&    target,
                MPI_Comm const& communicator,
                int             root = 0) {
  int stringSize, rank;
  MPI_Comm_rank(communicator, &rank);

  if (rank == root) {
    stringSize = target.size();
  }
  MPI_Bcast(&stringSize, 1, MPI_INT, 0, communicator);
  char* name = new char[stringSize + 1];

  if (rank == root) {
    target.copy(name, stringSize, 0);
  }
  name[stringSize] = '\0';
  MPI_Bcast(name, stringSize + 1, MPI_CHAR, 0, communicator);

  if (rank != root) {
    target = name;
  }
  delete[] name; name = NULL;
}

Configuration::
Configuration() {
  _filename = "";
}

Configuration::
Configuration(std::string const& filename) {
  _filename = filename;
}

void
Configuration::
setFileName(std::string const& filename) {
  _filename = filename;
}

void
Configuration::
loadParameters(Parameters&     parameters,
               MPI_Comm const& communicator) {
  tinyxml2::XMLDocument confFile;
  tinyxml2::XMLElement* node;
  tinyxml2::XMLElement* subNode;

  int rank;

  MPI_Comm_rank(communicator, &rank);

  if (rank == 0) {
    confFile.LoadFile(_filename.c_str());

    if (confFile.FirstChildElement() == NULL) {
      logError("Error parsing the configuration file");
    }

    node = confFile.FirstChildElement()->FirstChildElement("geometry");

    if (node == NULL) {
      logError("Error loading geometry properties");
    }

    readIntMandatory(parameters.geometry.sizeX, node, "sizeX");
    readIntMandatory(parameters.geometry.sizeY, node, "sizeY");
    readIntOptional(parameters.geometry.sizeZ, node, "sizeZ");

    if (parameters.geometry.sizeX < 2 || parameters.geometry.sizeY < 2 ||
        parameters.geometry.sizeZ < 0) {
      logError("Invalid size specified in configuration file");
    }

    parameters.geometry.dim = 0;

    if (node->QueryIntAttribute("dim", &(parameters.geometry.dim)) !=
        tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
      if (parameters.geometry.dim == 0) {
        if (parameters.geometry.sizeZ == 0) {
          parameters.geometry.sizeZ = 1;
          parameters.geometry.dim   = 2;
        } else {
          parameters.geometry.dim = 3;
        }
      }
    }

    if (parameters.geometry.dim == 3 && parameters.geometry.sizeZ == 1) {
      logError("Inconsistent data: 3D geometry specified with Z size zero");
    }

    readFloatMandatory(parameters.geometry.lengthX, node, "lengthX");
    readFloatMandatory(parameters.geometry.lengthY, node, "lengthY");
    readFloatMandatory(parameters.geometry.lengthZ, node, "lengthZ");

    std::string meshsizeType = "";
    subNode = node->FirstChildElement("mesh");
    readStringMandatory(meshsizeType, subNode);

    if (meshsizeType == "uniform") {
      parameters.geometry.meshsizeType = 0;
    } else if (meshsizeType == "stretched") {
      parameters.geometry.meshsizeType = 1;
      bool buffer = false;
        readBoolMandatory(buffer, node, "stretchX");
      parameters.geometry.stretchX = (int)buffer;
        readBoolMandatory(buffer, node, "stretchY");
      parameters.geometry.stretchY = (int)buffer;

      if (parameters.geometry.dim == 3) {
        readBoolMandatory(buffer, node, "stretchZ");
        parameters.geometry.stretchZ = (int)buffer;
      } else {
        parameters.geometry.stretchZ = false;
      }
    } else {
      logError("Unknown 'mesh'!");
    }

    _dim = parameters.geometry.dim;

    node = confFile.FirstChildElement()->FirstChildElement("timestep");

    if (node == NULL) {
      logError("Error loading timestep parameters");
    }

    readFloatOptional(parameters.timestep.dt,  node, "dt",    1);
    readFloatOptional(parameters.timestep.tau, node, "tau", 0.5);

    node = confFile.FirstChildElement()->FirstChildElement("flow");

    if (node == NULL) {
      logError("Error loading flow parameters");
    }

    readFloatMandatory(parameters.flow.Re, node, "Re");

    node = confFile.FirstChildElement()->FirstChildElement("solver");

    if (node == NULL) {
      logError("Error loading solver parameters");
    }

    readFloatMandatory(parameters.solver.gamma, node, "gamma");
    readIntOptional(parameters.solver.maxIterations, node, "maxIterations");

    node = confFile.FirstChildElement()->FirstChildElement("environment");

    if (node == NULL) {
      logError("Error loading environmental  parameters");
    }

    readFloatOptional(parameters.environment.gx, node, "gx");
    readFloatOptional(parameters.environment.gy, node, "gy");
    readFloatOptional(parameters.environment.gz, node, "gz");

    node = confFile.FirstChildElement()->FirstChildElement("simulation");

    if (node == NULL) {
      logError("Error loading simulation parameters");
    }

    readFloatMandatory(parameters.simulation.finalTime, node, "finalTime");

    subNode = node->FirstChildElement("type");

    if (subNode != NULL) {
      readStringMandatory(parameters.simulation.type, subNode);
    } else {
      logError("Missing type in simulation parameters");
    }

    subNode = node->FirstChildElement("scenario");

    if (subNode != NULL) {
      readStringMandatory(parameters.simulation.scenario, subNode);
    } else {
      logError("Missing scenario in simulation parameters");
    }

    node = confFile.FirstChildElement()->FirstChildElement("vtk");

    if (node == NULL) {
      logError("Error loading VTK parameters");
    }

    readFloatOptional(parameters.vtk.interval, node, "interval");
    readStringMandatory(parameters.vtk.prefix, node);

    node = confFile.FirstChildElement()->FirstChildElement("stdOut");

    if (node == NULL) {
      logError("Error loading StdOut parameters");
    }

    readFloatOptional(parameters.stdOut.interval, node, "interval", 1);

    node = confFile.FirstChildElement()->FirstChildElement("parallel");

    if (node == NULL) {
      logError("Error loading parallel parameters");
    }

    readIntOptional(parameters.parallel.numProcessors[0], node,
                    "numProcessorsX", 1);
    readIntOptional(parameters.parallel.numProcessors[1], node,
                    "numProcessorsY", 1);
    readIntOptional(parameters.parallel.numProcessors[2], node,
                    "numProcessorsZ", 1);

    parameters.parallel.leftNb   = MPI_PROC_NULL;
    parameters.parallel.rightNb  = MPI_PROC_NULL;
    parameters.parallel.bottomNb = MPI_PROC_NULL;
    parameters.parallel.topNb    = MPI_PROC_NULL;
    parameters.parallel.frontNb  = MPI_PROC_NULL;
    parameters.parallel.backNb   = MPI_PROC_NULL;

    parameters.parallel.localSize[0] = parameters.geometry.sizeX;
    parameters.parallel.localSize[1] = parameters.geometry.sizeY;
    parameters.parallel.localSize[2] = parameters.geometry.sizeZ;

    parameters.parallel.firstCorner[0] = 0;
    parameters.parallel.firstCorner[1] = 0;
    parameters.parallel.firstCorner[2] = 0;

    parameters.parallel.rank = rank;

    node = confFile.FirstChildElement()->FirstChildElement("walls");

    if (node == NULL) {
      logError("Error loading wall parameters");
    }

    tinyxml2::XMLElement* wall;
    wall = node->FirstChildElement("left");

    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorLeft, parameters.walls.scalarLeft);
    }

    wall = node->FirstChildElement("right");

    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorRight,
               parameters.walls.scalarRight);
    }

    wall = node->FirstChildElement("bottom");

    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorBottom,
               parameters.walls.scalarBottom);
    }

    wall = node->FirstChildElement("top");

    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorTop, parameters.walls.scalarTop);
    }

    wall = node->FirstChildElement("front");

    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorFront,
               parameters.walls.scalarFront);
    }

    wall = node->FirstChildElement("back");

    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorBack, parameters.walls.scalarBack);
    }

    if (parameters.simulation.scenario != "pressure-channel") {
      parameters.walls.scalarLeft = 0.0;
    }
    parameters.walls.scalarRight  = 0.0;
    parameters.walls.scalarBottom = 0.0;
    parameters.walls.scalarTop    = 0.0;
    parameters.walls.scalarFront  = 0.0;
    parameters.walls.scalarBack   = 0.0;

    parameters.bfStep.xRatio = -1.0;
    parameters.bfStep.yRatio = -1.0;
    node                     = confFile.FirstChildElement()->FirstChildElement(
      "backwardFacingStep");

    if (node != NULL) {
      readFloatMandatory(parameters.bfStep.xRatio, node, "xRatio");
      readFloatMandatory(parameters.bfStep.yRatio, node, "yRatio");
    }

    node = confFile.FirstChildElement()->FirstChildElement("turbulence_model");

    if (node != NULL) {
      readFloatOptional(parameters.blm.kappa,   node, "kappa");
      readFloatOptional(parameters.blm.delta99, node, "delta99");
      readStringMandatory(parameters.blm.modelType, node);
    }
  }

  MPI_Bcast(&(parameters.geometry.sizeX),        1, MPI_INT,    0,
            communicator);
  MPI_Bcast(&(parameters.geometry.sizeY),        1, MPI_INT,    0,
            communicator);
  MPI_Bcast(&(parameters.geometry.sizeZ),        1, MPI_INT,    0,
            communicator);

  MPI_Bcast(&(parameters.geometry.dim),          1, MPI_INT,    0,
            communicator);

  MPI_Bcast(&(parameters.geometry.meshsizeType), 1, MPI_INT,    0,
            communicator);
  MPI_Bcast(&(parameters.geometry.stretchX),     1, MPI_INT,    0,
            communicator);
  MPI_Bcast(&(parameters.geometry.stretchY),     1, MPI_INT,    0,
            communicator);
  MPI_Bcast(&(parameters.geometry.stretchZ),     1, MPI_INT,    0,
            communicator);
  MPI_Bcast(&(parameters.geometry.lengthX),      1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.geometry.lengthY),      1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.geometry.lengthZ),      1, MPI_DOUBLE, 0,
            communicator);

  MPI_Bcast(&(parameters.timestep.dt),           1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.timestep.tau),          1, MPI_DOUBLE, 0,
            communicator);

  MPI_Bcast(&(parameters.flow.Re),               1, MPI_DOUBLE, 0,
            communicator);

  MPI_Bcast(&(parameters.solver.gamma),          1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.solver.maxIterations),  1, MPI_DOUBLE, 0,
            communicator);

  MPI_Bcast(&(parameters.environment.gx),        1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.environment.gy),        1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.environment.gz),        1, MPI_DOUBLE, 0,
            communicator);

  MPI_Bcast(&(parameters.simulation.finalTime),  1, MPI_DOUBLE, 0,
            communicator);

  MPI_Bcast(&(parameters.vtk.interval),          1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.stdOut.interval),       1, MPI_INT,    0,
            communicator);

  broadcastString(parameters.vtk.prefix,          communicator);
  broadcastString(parameters.simulation.type,     communicator);
  broadcastString(parameters.simulation.scenario, communicator);

  MPI_Bcast(&(parameters.bfStep.xRatio),       1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.bfStep.yRatio),       1, MPI_DOUBLE, 0,
            communicator);

  MPI_Bcast(parameters.parallel.numProcessors, 3, MPI_INT,    0,
            communicator);

  MPI_Bcast(&(parameters.walls.scalarLeft),    1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.walls.scalarRight),   1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.walls.scalarBottom),  1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.walls.scalarTop),     1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.walls.scalarFront),   1, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(&(parameters.walls.scalarBack),    1, MPI_DOUBLE, 0,
            communicator);

  MPI_Bcast(parameters.walls.vectorLeft,       3, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(parameters.walls.vectorRight,      3, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(parameters.walls.vectorBottom,     3, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(parameters.walls.vectorTop,        3, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(parameters.walls.vectorFront,      3, MPI_DOUBLE, 0,
            communicator);
  MPI_Bcast(parameters.walls.vectorBack,       3, MPI_DOUBLE, 0,
            communicator);

  broadcastString(parameters.blm.modelType, communicator);
  MPI_Bcast(&(parameters.blm.delta99), 1, MPI_DOUBLE, 0, communicator);
  MPI_Bcast(&(parameters.blm.kappa),   1, MPI_DOUBLE, 0, communicator);
}
