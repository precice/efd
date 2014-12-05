#include "VTKStencil.h"

// TODO Correct performance issues for big matrices still required

void
VTKStencil::
writeVTKHeader(std::ostream& file) const {
  file << "# vtk DataFile Version 2.0" << std::endl
       << "I need something to put here" << std::endl
       << "ASCII" << std::endl << std::endl;
}

void
VTKStencil::
writePoints(std::ostream& file) const {
  // Number of points in every direction
  int px = _parameters.parallel.localSize[0] + 1;
  int py = _parameters.parallel.localSize[1] + 1;
  int pz = _parameters.geometry.dim == 2 ? 1 :
           _parameters.parallel.localSize[2] + 1;

  std::string grid;
  char        buffer[256];

  grid.reserve((file.precision() + 6) * px * py * pz * 3);

          sprintf(buffer,
          "DATASET STRUCTURED_GRID\nDIMENSIONS %d %d %d\nPOINTS %d float\n",
          px, py, pz, px * py * pz);
  grid.append(buffer);

  std::stringstream pointsStream; // ! Stream for the pressure data

  if (_parameters.geometry.dim == 3) {
    for (int k = 2; k < 2 + pz; k++) {
      for (int j = 2; j < 2 + py; j++) {
        for (int i = 2; i < 2 + px; i++) {
          // determine positions of grid points at lower/left/front corner of
          // the respective
          // grid cell (i,j,k) -> use meshsize-ptr
          sprintf(buffer, "%f %f %f\n", _parameters.meshsize->getPosX(i, j, k),
                  _parameters.meshsize->getPosY(i, j, k),
                  _parameters.meshsize->getPosZ(i, j, k));
          grid.append(buffer);
        }
      }
    }
  } else {
    for (int j = 2; j < 2 + py; j++) {
      for (int i = 2; i < 2 + px; i++) {
        // sprintf(buffer, "%f %f 0.0\n", _parameters.meshsize->getPosX(i, j),
        // _parameters.meshsize->getPosY(i, j));
        // grid.append(buffer);
        pointsStream << _parameters.meshsize->getPosX(i, j) << " "
                     << _parameters.meshsize->getPosY(i, j) << " "
                     << "0" << std::endl;
      }
    }
  }
  // grid.append("\n");

  file << grid;
  file << pointsStream.str() << std::endl;
}

VTKStencil::
VTKStencil(const Parameters& parameters)
  : FieldStencil<FlowField>(parameters), _prefix(parameters.vtk.prefix),
    _written(false) {}

void
VTKStencil::
apply(FlowField& flowField, int i, int j) {
  assertion(FieldStencil<FlowField>::_parameters.geometry.dim == 2);

  FLOAT pressure;
  FLOAT velocity[2] = { 0.0, 0.0 };

  if ((flowField.getFlags().getValue(i, j) & OBSTACLE_SELF) == 0) {
    flowField.getPressureAndVelocity(pressure, velocity, i, j);

    auto fgh = flowField.getFGH().getVector(i, j);
    //auto rhs = flowField.getRHS().getScalar(i, j);

    pressureStream << pressure << std::endl;

    velocityStream << velocity[0] << " "
                   << velocity[1] << " 0" << std::endl;

    fghStream << fgh[0] << " "
              << fgh[1] << " "
              << "0" << std::endl;

    //rhsStream << rhs << std::endl;
  } else {
    pressureStream << "0.0" << std::endl;
    velocityStream << "0.0 0.0 0.0" << std::endl;
  }
}

void
VTKStencil::
apply(FlowField& flowField, int i, int j, int k) {
  assertion(FieldStencil<FlowField>::_parameters.geometry.dim == 3);

  FLOAT pressure;
  FLOAT velocity[3] = { 0.0, 0.0, 0.0 };

  if ((flowField.getFlags().getValue(i, j, k) & OBSTACLE_SELF) == 0) {
    flowField.getPressureAndVelocity(pressure, velocity, i, j, k);

    pressureStream << pressure << std::endl;

    velocityStream << velocity[0] << " "
                   << velocity[1] << " "
                   << velocity[2] << std::endl;
  } else {
    pressureStream << "0.0" << std::endl;
    velocityStream << "0.0 0.0 0.0" << std::endl;
  }
}

void
VTKStencil::
openFile(const FlowField& flowField, int timeStep) {
  // This relates to openning the file with the correct name
  _written = false;
  std::stringstream namestream;
  std::string       name;
  namestream.precision(4);
  namestream << _prefix << "." << _parameters.parallel.rank << "." <<
    timeStep << ".vtk";
  name = namestream.str();
  _ofile.open(name.c_str());
  namestream.str("");

  writeVTKHeader(_ofile);

  writePoints(_ofile);
}

void
VTKStencil::
write(FlowField& flowField, int timeStep) {
  openFile(flowField, timeStep);

  if (FieldStencil<FlowField>::_parameters.geometry.dim == 2) {
    // Write pressure
    _ofile << "CELL_DATA " << flowField.getNx() * flowField.getNy() << std::endl
           << "SCALARS pressure float 1" << std::endl
           << "LOOKUP_TABLE default" << std::endl;
    _ofile << fghStream.str() << std::endl;
    fghStream.str("");
    //_ofile << rhsStream.str() << std::endl;
    rhsStream.str("");
    // Write velocity
    // _ofile << "CELL_DATA " << flowField.getNx() * flowField.getNy() <<
    // std::endl
    // _ofile << "VECTORS velocity float" << std::endl;
    _ofile << velocityStream.str() << std::endl;
    velocityStream.str("");
    _ofile << pressureStream.str() << std::endl;
    pressureStream.str("");
  }

  if (FieldStencil<FlowField>::_parameters.geometry.dim == 3) {
    // Write pressure
    _ofile << "CELL_DATA " << flowField.getNx() * flowField.getNy() *
      flowField.getNz()
           << std::endl
           << "SCALARS pressure float 1" << std::endl
           << "LOOKUP_TABLE default" << std::endl;
    _ofile << pressureStream.str() << std::endl;
    pressureStream.str("");

    // Write velocity
    // _ofile << "CELL_DATA " << flowField.getNx() * flowField.getNy() <<
    // std::endl
    _ofile << "VECTORS velocity float" << std::endl;
    _ofile << velocityStream.str() << std::endl;
    velocityStream.str("");
  }

  _written = true;
  closeFile();
}

void
VTKStencil::
closeFile() {
  _ofile.close();
}
