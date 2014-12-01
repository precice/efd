#ifndef FsiSimulation_mpigenerics_hpp
#define FsiSimulation_mpigenerics_hpp

#include <mpi.h>

namespace FsiSimulation {
template <typename Scalar>
MPI_Datatype
getMpiScalarType() { return MPI_FLOAT; }

template <>
MPI_Datatype
getMpiScalarType<float
                 >() { return MPI_FLOAT; }

template <>
MPI_Datatype
getMpiScalarType<double
                 >() { return MPI_DOUBLE; }
}

#endif
