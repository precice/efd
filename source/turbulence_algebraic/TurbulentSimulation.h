#ifndef _TURBULENTSIMULATION_H_
#define _TURBULENTSIMULATION_H_

#include "../Simulation.h"
#include "TurbulentFlowField.h"
#include "stencils/ComputeLocalViscosityStencil.h"
#include "stencils/TurbFGHStencil.h"
#include "stencils/VTKStencil4Turbulence.h"
#include "stencils/MaximumTurbulentViscosityStencil.h"
#include "parallelManagers/TurbulentPetscParallelManager.h"
#include "TurbulenceViscosityFactory.h"

/** inherits from the Simulation class and performs one time step, incorporating the
 *  functionality for turbulence simulations.
 *
 *  @author Philipp Neumann
 */
class TurbulentSimulation: public Simulation {
private:
	TurbulenceViscosityFactory _stencilFactory;
	TurbulentFlowField &_turbulentFlowField; //! New data structure for turbulence simulation
	ComputeLocalViscosityStencil *_computeLocalViscosityStencil; //! stencil for local viscosity evaluations
	FieldIterator<TurbulentFlowField> _computeLocalViscosityIterator; //! iterator for the local viscosity evaluation stencil
	MaximumTurbulentViscosityStencil _maxTurbViscStencil; //! stencil for computation of minimal turb. viscosity
	FieldIterator<TurbulentFlowField> _maxTurbViscIterator; //! iterator for the stencil
	TurbFGHStencil _turbulentFGHStencil;
	FieldIterator<TurbulentFlowField> _turbulentFGHIterator;
	TurbulentPetscParallelManager _turbulentParallelManager;

        /** initialises the fields for turbulence simulation (viscosity, distance from wall) */
        void initializeTurbulenceFlowfield();



public:
	TurbulentSimulation(Parameters &parameters, TurbulentFlowField &flowField);
	virtual ~TurbulentSimulation();

	virtual void solveTimestep();

	virtual void initializeFlowField(){
	  // first: initialize step-geometry and other scenario-dependent things
	  Simulation::initializeFlowField();
	  // second: set up fields for turbulence simulation
	  initializeTurbulenceFlowfield();
	}

    virtual void plotVTK(int timeStep){
      // use the new plotter stencil to also visualise the viscosity field
      VTKStencil4Turbulence vtkStencil( _parameters );
      FieldIterator<TurbulentFlowField> vtkIterator( _turbulentFlowField, _parameters,vtkStencil, 1, 0 );

      vtkIterator.iterate();
      vtkStencil.write( _turbulentFlowField,timeStep );
    }

  protected:
    virtual void setTimeStep(){
      FLOAT localMin, globalMin;
      FLOAT factor = 1.0/(_parameters.meshsize->getDxMin() * _parameters.meshsize->getDxMin()) +
                     1.0/(_parameters.meshsize->getDyMin() * _parameters.meshsize->getDyMin());

      // compute the maximum velocity first
      _maxUStencil.reset();
      _maxUFieldIterator.iterate();
      _maxUBoundaryIterator.iterate();
      // compute the maximum turb. viscosity
      _maxTurbViscStencil.reset();
      _maxTurbViscIterator.iterate();

      if (_parameters.geometry.dim == 3) {
          factor += 1.0/(_parameters.meshsize->getDzMin() * _parameters.meshsize->getDzMin());
          _parameters.timestep.dt = 1.0 / _maxUStencil.getMaxValues()[2];
      } else {
          _parameters.timestep.dt = 1.0 / _maxUStencil.getMaxValues()[0];
      }

      // incorporate minimum turbulent viscosity into choice of time step
      localMin = std::min(_parameters.timestep.dt,
                          std::min(std::min( 1/(_maxTurbViscStencil.getMaximumTurbulentViscosity()+1/_parameters.flow.Re)/(2*factor),
                          1 / _maxUStencil.getMaxValues()[0]),
                          1 / _maxUStencil.getMaxValues()[1]));

      // Here, we select the type of operation before compiling. This allows to use the correct
      // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
      // machines.

      globalMin = MY_FLOAT_MAX;
      MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

      _parameters.timestep.dt = globalMin;
      _parameters.timestep.dt *= _parameters.timestep.tau;
    }

};

#endif // _TURBULENTSIMULATION_H_
