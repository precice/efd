#ifndef _PETSC_PARALLEL_MANAGER_H_
#define _PETSC_PARALLEL_MANAGER_H_

#include <petscdmda.h>
#include <petscksp.h>

#include "../Iterators.h"
#include "BufferStencils.h"
#include "FlowField.h"
#include "Parameters.h"
#include "Solvers/PetscSolver.hpp"

/** Class dedicated to the transfer of boundary information between subdomains
 */
class PetscParallelManager {
protected:
  FlowField&        _flowField;
  const Parameters& _parameters;

  // Location of the elements of the subdomain
  PetscInt _firstX, _lengthX, _firstY, _lengthY, _firstZ, _lengthZ;

  // Input and output buffers
  FLOAT* _leftBufferIn, * _rightBufferIn,
  * _bottomBufferIn, * _topBufferIn,
  * _frontBufferIn, * _backBufferIn;

  FLOAT* _leftBufferOut, * _rightBufferOut,
  * _bottomBufferOut, * _topBufferOut,
  * _frontBufferOut, * _backBufferOut;

  int _bufferSize[3]; // ! Sizes of the buffers in each direction.
                      // Components are x, y, z

  int _pressureSize[3];         // ! Size for pressure transfer. Smaller than
                                // required for velocity
  MPI_Status _mpiStatus; // ! For the MPI functions

private:
  // ! @brief Functions to fill and read the buffers
  // ! @{
  VelocityBufferFillStencil _fillVelocityStencil;
  VelocityBufferReadStencil _readVelocityStencil;

  PressureBufferFillStencil _fillPressureStencil;
  PressureBufferReadStencil _readPressureStencil;

  ParallelBoundaryIterator<FlowField> _fillVelocityIterator;
  ParallelBoundaryIterator<FlowField> _readVelocityIterator;

  ParallelBoundaryIterator<FlowField> _fillPressureIterator;
  ParallelBoundaryIterator<FlowField> _readPressureIterator;
  // !@}

public:
  /* Constructor
   * @param parameters An instance of the parameters
   */
  PetscParallelManager(FlowField& flowField, const Parameters& parameters);

  /** Destructor */
  virtual
  ~PetscParallelManager();

  /** Communicates boundary values for pressure */
  void
  communicatePressure();

  /** Communicates boundary values of velocity */
  void
  communicateVelocity();
};

#endif
