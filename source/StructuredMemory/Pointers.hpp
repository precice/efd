#ifndef FsiSimulation_StructuredMemory_Pointers_hpp
#define FsiSimulation_StructuredMemory_Pointers_hpp

#include <Eigen/Core>

namespace FsiSimulation {
namespace StructuredMemory {
template <typename TCell, int D>
struct Pointers {};

template <typename TCell>
struct Pointers<TCell, 1> {
  typedef TCell*                   Type;
  typedef Eigen::Matrix<int, 1, 1> VectorDi;

  Pointers() : _pointers(0) {}

  Pointers(Pointers const& other) : _size(other._size),
                                    _pointers(other._pointers) {}

  ~Pointers() {
    release();
  }

  Pointers&
  operator=(Pointers const& other) {
    _size     = other._size;
    _pointers = other._pointers;

    return *this;
  }

  TCell*
  getCell(VectorDi const& index) {
    return &_pointers[index(0)];
  }

  void
  allocate(VectorDi const& size) {
    release();
    _size     = size;
    _pointers = new TCell[_size(0)];
  }

  void
  release() {
    if (_pointers) {
      delete[] _pointers;
    }
  }

  static inline TCell&
  dereference(Type _pointers, VectorDi const& index) {
    return _pointers[index(0)];
  }

private:
  VectorDi _size;
  Type     _pointers;
};

template <typename TCell>
struct Pointers<TCell, 2> {
  typedef TCell**                  Type;
  typedef Eigen::Matrix<int, 2, 1> VectorDi;

  Pointers() : _data(0),
               _pointers(0) {}

  Pointers(Pointers const& other) : _size(other._size),
                                    _data(other._data),
                                    _pointers(other._pointers) {}

  ~Pointers() {
    release();
  }

  Pointers&
  operator=(Pointers const& other) {
    _size     = other._size;
    _data     = other._data;
    _pointers = other._pointers;

    return *this;
  }

  TCell*
  getCell(VectorDi const& index) {
    return &_pointers[index(1)][index(0)];
  }

  void
  allocate(VectorDi const& size) {
    release();
    _size     = size;
    _data     = new TCell[_size.prod()];
    _pointers = new TCell*[_size(1)];

    for (int j = 0; j < _size(1); ++j) {
      _pointers[j] = &_data[_size(0) * j];
    }
  }

  void
  release() {
    if (_data) {
      delete[] _data;
    }

    if (_pointers) {
      delete[] _pointers;
    }
  }

  static inline TCell&
  dereference(Type pointers, VectorDi const& index) {
    return pointers[index(1)][index(0)];
  }

private:
  VectorDi _size;
  TCell*   _data;
  Type     _pointers;
};

template <typename TCell>
struct Pointers<TCell, 3> {
  typedef TCell***                 Type;
  typedef Eigen::Matrix<int, 3, 1> VectorDi;

  Pointers() : _data(0),
               _pointers(0) {}

  Pointers(Pointers const& other) : _size(other._size),
                                    _data(other._data),
                                    _pointers(other._pointers) {}

  ~Pointers() {
    release();
  }

  Pointers&
  operator=(Pointers const& other) {
    _size     = other._size;
    _data     = other._data;
    _pointers = other._pointers;

    return *this;
  }

  TCell*
  getCell(VectorDi const& index) {
    return &_pointers[index(2)][index(1)][index(0)];
  }

  void
  allocate(VectorDi const& size) {
    release();
    _size     = size;
    _data     = new TCell[_size.prod()];
    _pointers = new TCell * *[_size(2)];

    for (int k = 0; k < _size(2); ++k) {
      _pointers[k] = new TCell*[_size(1)];

      for (int j = 0; j < _size(1); ++j) {
        _pointers[j][k] =
          &_data[_size(1) * _size(0) * k + _size(0) * j];
      }
    }
  }

  void
  release() {
    if (_data) {
      delete[] _data;
    }

    if (_pointers) {
      for (int k = 0; k < _size(2); ++k) {
        delete[] _pointers[k];
      }
      delete[] _pointers;
    }
  }

  static inline TCell&
  dereference(Type pointers, VectorDi const& index) {
    return pointers[index(2)][index(1)][index(0)];
  }

private:
  VectorDi _size;
  TCell*   _data;
  Type     _pointers;
};
}
}
#endif
