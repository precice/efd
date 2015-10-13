#ifndef FsiSimulation_FluidSimulation_Private_mpigenerics_hpp
#define FsiSimulation_FluidSimulation_Private_mpigenerics_hpp

#include <Uni/ExecutionControl/exception>

#include <mpi.h>

namespace FsiSimulation {
namespace FluidSimulation {
namespace Private {
template <typename TScalar>
inline MPI_Datatype
getMpiScalarType() {
  throwException("This type is not supported for getMpiScalarType()");

  return MPI_FLOAT;
}

template <>
inline MPI_Datatype
getMpiScalarType<int>() {
  return MPI_INT;
}

template <>
inline MPI_Datatype
getMpiScalarType<unsigned>() {
  return MPI_UNSIGNED;
}

template <>
inline MPI_Datatype
getMpiScalarType<float>() {
  return MPI_FLOAT;
}

template <>
inline MPI_Datatype
getMpiScalarType<double>() {
  return MPI_DOUBLE;
}

template <>
inline MPI_Datatype
getMpiScalarType<long double>() {
  return MPI_LONG_DOUBLE;
}

template <typename TScalar>
inline void
mpi_all_reduce(void*        sendbuf,
               TScalar*        recvbuf,
               int const&      count,
               MPI_Op const&   op,
               MPI_Comm const& comm) {
  MPI_Allreduce(sendbuf, recvbuf, count, getMpiScalarType<TScalar>(), op, comm);
}
}
}
}

#endif
