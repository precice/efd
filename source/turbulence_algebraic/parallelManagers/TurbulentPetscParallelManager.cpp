#include "TurbulentPetscParallelManager.h"
//

TurbulentPetscParallelManager::
TurbulentPetscParallelManager(TurbulentFlowField& flowField, const
                              Parameters& parameters)
  : PetscParallelManager(flowField, parameters),
    _fillViscosityStencil(parameters, _leftBufferOut, _rightBufferOut,
                          _bottomBufferOut,
                          _topBufferOut, _frontBufferOut, _backBufferOut),
    _readViscosityStencil(parameters, _leftBufferIn, _rightBufferIn,
                          _bottomBufferIn,
                          _topBufferIn, _frontBufferIn, _backBufferIn),
    _fillViscosityIterator(flowField, parameters, _fillViscosityStencil, 1, 0),
    _readViscosityIterator(flowField, parameters, _readViscosityStencil, 1,
                           0) {}

TurbulentPetscParallelManager::~TurbulentPetscParallelManager() {}

void
TurbulentPetscParallelManager::
communicateViscosity() {
  _fillViscosityIterator.iterate();
    MPI_Sendrecv(PetscParallelManager::_rightBufferOut,
               PetscParallelManager::_pressureSize[0],
               MY_MPI_FLOAT, PetscParallelManager::_parameters.parallel.rightNb,
               0,
               PetscParallelManager::_leftBufferIn,
               PetscParallelManager::_pressureSize[0],
               MY_MPI_FLOAT, PetscParallelManager::_parameters.parallel.leftNb,
               0,
               PETSC_COMM_WORLD, &(PetscParallelManager::_mpiStatus));
    MPI_Sendrecv(PetscParallelManager::_leftBufferOut,
               PetscParallelManager::_pressureSize[0],
               MY_MPI_FLOAT, PetscParallelManager::_parameters.parallel.leftNb,
               1,
               PetscParallelManager::_rightBufferIn,
               PetscParallelManager::_pressureSize[0],
               MY_MPI_FLOAT, PetscParallelManager::_parameters.parallel.rightNb,
               1,
               PETSC_COMM_WORLD, &(PetscParallelManager::_mpiStatus));
    MPI_Sendrecv(PetscParallelManager::_topBufferOut,
               PetscParallelManager::_pressureSize[1],
               MY_MPI_FLOAT, PetscParallelManager::_parameters.parallel.topNb,
               2,
               PetscParallelManager::_bottomBufferIn,
               PetscParallelManager::_pressureSize[1],
               MY_MPI_FLOAT,
               PetscParallelManager::_parameters.parallel.bottomNb, 2,
               PETSC_COMM_WORLD, &(PetscParallelManager::_mpiStatus));
    MPI_Sendrecv(PetscParallelManager::_bottomBufferOut,
               PetscParallelManager::_pressureSize[1],
               MY_MPI_FLOAT,
               PetscParallelManager::_parameters.parallel.bottomNb, 3,
               PetscParallelManager::_topBufferIn,
               PetscParallelManager::_pressureSize[1],
               MY_MPI_FLOAT, PetscParallelManager::_parameters.parallel.topNb,
               3,
               PETSC_COMM_WORLD, &(PetscParallelManager::_mpiStatus));

  if (PetscParallelManager::_parameters.geometry.dim == 3) {
    MPI_Sendrecv(PetscParallelManager::_backBufferOut,
                 PetscParallelManager::_pressureSize[2],
                 MY_MPI_FLOAT,
                 PetscParallelManager::_parameters.parallel.backNb, 4,
                 PetscParallelManager::_frontBufferIn,
                 PetscParallelManager::_pressureSize[2],
                 MY_MPI_FLOAT,
                 PetscParallelManager::_parameters.parallel.frontNb, 4,
                 PETSC_COMM_WORLD, &(PetscParallelManager::_mpiStatus));
    MPI_Sendrecv(PetscParallelManager::_frontBufferOut,
                 PetscParallelManager::_pressureSize[2],
                 MY_MPI_FLOAT,
                 PetscParallelManager::_parameters.parallel.frontNb, 5,
                 PetscParallelManager::_backBufferIn,
                 PetscParallelManager::_pressureSize[2],
                 MY_MPI_FLOAT,
                 PetscParallelManager::_parameters.parallel.backNb, 5,
                 PETSC_COMM_WORLD, &(PetscParallelManager::_mpiStatus));
  }
  _readViscosityIterator.iterate();
}
