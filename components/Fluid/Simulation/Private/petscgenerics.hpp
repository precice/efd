#ifndef FsiSimulation_FluidSimulation_Private_petscgenerics_hpp
#define FsiSimulation_FluidSimulation_Private_petscgenerics_hpp

#include <Eigen/Core>

#include <memory>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

namespace FsiSimulation {
namespace FluidSimulation {
#if ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 5))
template <int TD>
using DMBoundaryTypeVector = Eigen::Matrix<DMBoundaryType, TD, 1>;
#else
template <int TD>
using DMBoundaryTypeVector = Eigen::Matrix<DMDABoundaryType, TD, 1>;
#endif

template <int TD>
using VectorDi = Eigen::Matrix<int, TD, 1>;

typedef std::unique_ptr<PetscInt const[]> UniqueConstPetscIntArray;
template <int TD>
using VectorDConstPetscIntPointer
        = Eigen::Matrix<UniqueConstPetscIntArray, TD, 1>;

template <int TD>
inline DMBoundaryTypeVector<TD>
createDMBoundaries() {
  DMBoundaryTypeVector<TD> result;

  for (int i = 0; i < TD; ++i) {
#if ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 5))
    result(i) = DM_BOUNDARY_NONE;
#else
    result(i) = DMDA_BOUNDARY_NONE;
#endif
  }

  return result;
}

template <int TD>
inline void
DMCreate(MPI_Comm const&,
         DMBoundaryTypeVector<TD> const&,
         DMDAStencilType const&,
         VectorDi<TD> const&,
         VectorDi<TD> const&,
         int const&,
         int const&,
         VectorDConstPetscIntPointer<TD> const&,
         DM*) {}

template <>
inline void
DMCreate<2
         >(MPI_Comm const&                       comm,
           DMBoundaryTypeVector<2> const&        boundaryTypes,
           DMDAStencilType const&                stencilType,
           VectorDi<2> const&                    globalSize,
           VectorDi<2> const&                    processorSize,
           int const&                            dof,
           int const&                            stencilWidth,
           VectorDConstPetscIntPointer<2> const& localSizes,
           DM*                                   da) {
  DMDACreate2d(comm,
               boundaryTypes(0),
               boundaryTypes(1),
               stencilType,
               globalSize(0),
               globalSize(1),
               processorSize(0),
               processorSize(1),
               dof,
               stencilWidth,
               localSizes(0).get(),
               localSizes(1).get(),
               da);
}

template <>
inline void
DMCreate<3
         >(MPI_Comm const&                       comm,
           DMBoundaryTypeVector<3> const&        boundaryTypes,
           DMDAStencilType const&                stencilType,
           VectorDi<3> const&                    globalSize,
           VectorDi<3> const&                    processorSize,
           int const&                            dof,
           int const&                            stencilWidth,
           VectorDConstPetscIntPointer<3> const& localSizes,
           DM*                                   da) {
  DMDACreate3d(comm,
               boundaryTypes(0),
               boundaryTypes(1),
               boundaryTypes(2),
               stencilType,
               globalSize(0),
               globalSize(1),
               globalSize(2),
               processorSize(0),
               processorSize(1),
               processorSize(2),
               dof,
               stencilWidth,
               localSizes(0).get(),
               localSizes(1).get(),
               localSizes(2).get(),
               da);
}
}
}
#endif
