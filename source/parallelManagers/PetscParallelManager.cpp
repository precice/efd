#include "PetscParallelManager.h"
#include <iostream>
PetscParallelManager::
PetscParallelManager(FlowField& flowField, const Parameters& parameters)
  : _flowField(flowField), _parameters(parameters),

    _fillVelocityStencil(parameters, _leftBufferOut, _rightBufferOut,
                         _bottomBufferOut,
                         _topBufferOut, _frontBufferOut, _backBufferOut),
    _readVelocityStencil(parameters, _leftBufferIn, _rightBufferIn,
                         _bottomBufferIn,
                         _topBufferIn, _frontBufferIn, _backBufferIn),
    _fillPressureStencil(parameters, _leftBufferOut, _rightBufferOut,
                         _bottomBufferOut,
                         _topBufferOut, _frontBufferOut, _backBufferOut),
    _readPressureStencil(parameters, _leftBufferIn, _rightBufferIn,
                         _bottomBufferIn,
                         _topBufferIn, _frontBufferIn, _backBufferIn),

    _fillVelocityIterator(_flowField, parameters, _fillVelocityStencil, 1, 0),
    _readVelocityIterator(_flowField, parameters, _readVelocityStencil, 1, 0),
    _fillPressureIterator(_flowField, parameters, _fillPressureStencil, 1, 0),
    _readPressureIterator(_flowField, parameters, _readPressureStencil, 1, 0) {
  // Allocate buffers and set number of pressure values
  if (_parameters.geometry.dim == 2) {
    _bufferSize[0] = 2 * (flowField.getNy() + 3); // Allocate space for all
                                                  // values
    _bufferSize[1] = 2 * (flowField.getNx() + 3);

    _pressureSize[0] = (flowField.getNy() + 3);
    _pressureSize[1] = (flowField.getNx() + 3);

    _leftBufferIn  = new FLOAT[_bufferSize[0]];
    _leftBufferOut = new FLOAT[_bufferSize[0]];

    _rightBufferIn  = new FLOAT[_bufferSize[0]];
    _rightBufferOut = new FLOAT[_bufferSize[0]];

    _bottomBufferIn  = new FLOAT[_bufferSize[1]];
    _bottomBufferOut = new FLOAT[_bufferSize[1]];

    _topBufferIn  = new FLOAT[_bufferSize[1]];
    _topBufferOut = new FLOAT[_bufferSize[1]];
  } else if (_parameters.geometry.dim == 3) {
    _bufferSize[0] = 3 * ((flowField.getNy() + 3) * (flowField.getNz() + 3));
    _bufferSize[1] = 3 * ((flowField.getNx() + 3) * (flowField.getNz() + 3));
    _bufferSize[2] = 3 * ((flowField.getNx() + 3) * (flowField.getNy() + 3));

    _pressureSize[0] = (flowField.getNy() + 3) * (flowField.getNz() + 3);
    _pressureSize[1] = (flowField.getNx() + 3) * (flowField.getNz() + 3);
    _pressureSize[2] = (flowField.getNx() + 3) * (flowField.getNy() + 3);

    _leftBufferIn  = new FLOAT[_bufferSize[0]];
    _leftBufferOut = new FLOAT[_bufferSize[0]];

    _rightBufferIn  = new FLOAT[_bufferSize[0]];
    _rightBufferOut = new FLOAT[_bufferSize[0]];

    _bottomBufferIn  = new FLOAT[_bufferSize[1]];
    _bottomBufferOut = new FLOAT[_bufferSize[1]];

    _topBufferIn  = new FLOAT[_bufferSize[1]];
    _topBufferOut = new FLOAT[_bufferSize[1]];

    _frontBufferIn  = new FLOAT[_bufferSize[2]];
    _frontBufferOut = new FLOAT[_bufferSize[2]];

    _backBufferIn  = new FLOAT[_bufferSize[2]];
    _backBufferOut = new FLOAT[_bufferSize[2]];
  }
}

PetscParallelManager::~PetscParallelManager() {
  delete[] _leftBufferIn;
  delete[] _rightBufferIn;
  delete[] _bottomBufferIn;
  delete[] _topBufferIn;
  delete[] _leftBufferOut;
  delete[] _rightBufferOut;
  delete[] _bottomBufferOut;
  delete[] _topBufferOut;

  if (_parameters.geometry.dim == 3) {
    delete[] _frontBufferIn;
    delete[] _frontBufferOut;
    delete[] _backBufferIn;
    delete[] _backBufferOut;
  }
}

void
PetscParallelManager::
communicatePressure() {
  _fillPressureIterator.iterate();
    MPI_Sendrecv(_rightBufferOut, _pressureSize[0], MY_MPI_FLOAT,
               _parameters.parallel.rightNb,  0,
               _leftBufferIn, _pressureSize[0], MY_MPI_FLOAT,
               _parameters.parallel.leftNb, 0,
               PETSC_COMM_WORLD, &_mpiStatus);
    MPI_Sendrecv(_leftBufferOut, _pressureSize[0], MY_MPI_FLOAT,
               _parameters.parallel.leftNb,   1,
               _rightBufferIn, _pressureSize[0], MY_MPI_FLOAT,
               _parameters.parallel.rightNb, 1,
               PETSC_COMM_WORLD, &_mpiStatus);
    MPI_Sendrecv(_topBufferOut, _pressureSize[1], MY_MPI_FLOAT,
               _parameters.parallel.topNb,    2,
               _bottomBufferIn, _pressureSize[1], MY_MPI_FLOAT,
               _parameters.parallel.bottomNb, 2,
               PETSC_COMM_WORLD, &_mpiStatus);
    MPI_Sendrecv(_bottomBufferOut, _pressureSize[1], MY_MPI_FLOAT,
               _parameters.parallel.bottomNb, 3,
               _topBufferIn, _pressureSize[1], MY_MPI_FLOAT,
               _parameters.parallel.topNb, 3,
               PETSC_COMM_WORLD, &_mpiStatus);

  if (_parameters.geometry.dim == 3) {
    MPI_Sendrecv(_backBufferOut, _pressureSize[2], MY_MPI_FLOAT,
                 _parameters.parallel.backNb,  4,
                 _frontBufferIn, _pressureSize[2], MY_MPI_FLOAT,
                 _parameters.parallel.frontNb, 4,
                 PETSC_COMM_WORLD, &_mpiStatus);
    MPI_Sendrecv(_frontBufferOut, _pressureSize[2], MY_MPI_FLOAT,
                 _parameters.parallel.frontNb, 5,
                 _backBufferIn, _pressureSize[2], MY_MPI_FLOAT,
                 _parameters.parallel.backNb, 5,
                 PETSC_COMM_WORLD, &_mpiStatus);
  }
  _readPressureIterator.iterate();
}

void
PetscParallelManager::
communicateVelocity() {
  _fillVelocityIterator.iterate();
    MPI_Sendrecv(_rightBufferOut, _bufferSize[0], MY_MPI_FLOAT,
               _parameters.parallel.rightNb,  0,
               _leftBufferIn, _bufferSize[0], MY_MPI_FLOAT,
               _parameters.parallel.leftNb, 0,
               PETSC_COMM_WORLD, &_mpiStatus);
    MPI_Sendrecv(_leftBufferOut, _bufferSize[0], MY_MPI_FLOAT,
               _parameters.parallel.leftNb,   1,
               _rightBufferIn, _bufferSize[0], MY_MPI_FLOAT,
               _parameters.parallel.rightNb, 1,
               PETSC_COMM_WORLD, &_mpiStatus);
    MPI_Sendrecv(_topBufferOut, _bufferSize[1], MY_MPI_FLOAT,
               _parameters.parallel.topNb,    2,
               _bottomBufferIn, _bufferSize[1], MY_MPI_FLOAT,
               _parameters.parallel.bottomNb, 2,
               PETSC_COMM_WORLD, &_mpiStatus);
    MPI_Sendrecv(_bottomBufferOut, _bufferSize[1], MY_MPI_FLOAT,
               _parameters.parallel.bottomNb, 3,
               _topBufferIn, _bufferSize[1], MY_MPI_FLOAT,
               _parameters.parallel.topNb, 3,
               PETSC_COMM_WORLD, &_mpiStatus);

  if (_parameters.geometry.dim == 3) {
    MPI_Sendrecv(_backBufferOut, _bufferSize[2], MY_MPI_FLOAT,
                 _parameters.parallel.backNb,  4,
                 _frontBufferIn, _bufferSize[2], MY_MPI_FLOAT,
                 _parameters.parallel.frontNb, 4,
                 PETSC_COMM_WORLD, &_mpiStatus);
    MPI_Sendrecv(_frontBufferOut, _bufferSize[2], MY_MPI_FLOAT,
                 _parameters.parallel.frontNb, 5,
                 _backBufferIn, _bufferSize[2], MY_MPI_FLOAT,
                 _parameters.parallel.backNb, 5,
                 PETSC_COMM_WORLD, &_mpiStatus);
  }
  _readVelocityIterator.iterate();
}
