#ifndef _VTKSTENCIL4TURBULENCE_H_
#define _VTKSTENCIL4TURBULENCE_H_

#include "Definitions.h"
#include "Parameters.h"
#include "Stencil.h"
#include "TurbulentFlowField.h"
#include <fstream>
#include <sstream>
#include <string>

class VTKStencil4Turbulence : public FieldStencil<TurbulentFlowField> {
private:
  void
  writeVTKHeader(std::ostream& file) const;
  void
  writePoints(std::ostream& file) const;

  std::string   _prefix;    // ! Prefix to be attached to the vtk files
  std::ofstream _ofile;     // ! Output file stream
  bool          _written;   // ! Whether the file has already been written

  std::stringstream pressureStream;     // ! Stream for the pressure data
  std::stringstream velocityStream;     // ! Stream for the velocity data
  std::stringstream viscosityStream;     // ! Stream for the velocity data

  /** Open output file
   * Opens the output file and prepares for writing.
   * It also writes the header of the file and the points og the grid.
   * @param flowField State of the flow field
   * @param time Current simulation time in double format
   */
  void
  openFile(const TurbulentFlowField& flowField, int timeStep);

  /** Finish writing. Must be called once the file has been written.
   *
   * Stores all the streams and closes the file.
   */
  void
  closeFile();

public:
  /** Constructor
   *
   * @param prefix String with the prefix of the name of the VTK files
   */
  VTKStencil4Turbulence(const Parameters& parameters);
  ~VTKStencil4Turbulence() {}

  /** 2D operation for one position
   *
   * @param flowField State of the flow field
   * @param i Position in the x direction
   * @param j Position in the y direction
   */
  virtual void
  apply(TurbulentFlowField& flowField, int i, int j);

  /** 3D operation for one position
   *
   * @param flowField State of the flow field
   * @param i Position in the x direction
   * @param j Position in the y direction
   * @param k Position in the z direction
   */
  virtual void
  apply(TurbulentFlowField& flowField, int i, int j, int k);

  /** Writes the information to the file
   * @param flowField Flow field to be written
   */
  virtual void
  write(TurbulentFlowField& flowField, int timeStep);
};

#endif // _VTKSTENCIL4TURBULENCE_H_
