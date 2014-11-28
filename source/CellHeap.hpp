#ifndef FsiSimulation_CellHeap_hpp
#define FsiSimulation_CellHeap_hpp

#include <Eigen/Core>

namespace FsiSimulation {
namespace Private {
template <typename Cell, int D>
struct Pointers {};

template <typename Cell>
struct Pointers<Cell, 1> {
  typedef Eigen::Matrix<int, 1, 1> VectorDi;

  Pointers() : pointers(0) {}

  Pointers(Pointers const& other) : size(other.size),
                                    pointers(other.pointers) {}

  ~Pointers() {
    release();
  }

  Pointers&
  operator=(Pointers const& other) {
    size     = other.size;
    pointers = other.pointers;

    return *this;
  }

  Cell*
  getCell(VectorDi const& index) {
    return &pointers[index(0)];
  }

  void
  allocate(VectorDi const& size_) {
    release();
    size     = size_;
    pointers = new Cell[size(0)];
  }

  void
  release() {
    if (pointers) {
      delete pointers;
    }
  }

  VectorDi size;
  Cell*    pointers;
};
template <typename Cell>
struct Pointers<Cell, 2> {
  typedef Eigen::Matrix<int, 2, 1> VectorDi;

  Pointers() : _data(0),
               pointers(0) {}

  Pointers(Pointers const& other) : size(other.size),
                                    _data(other._data),
                                    pointers(other.pointers) {}

  ~Pointers() {
    release();
  }

  Pointers&
  operator=(Pointers const& other) {
    size     = other.size;
    _data    = other._data;
    pointers = other.pointers;

    return *this;
  }

  Cell*
  getCell(VectorDi const& index) {
    return &pointers[index(1)][index(0)];
  }

  void
  allocate(VectorDi const& size_) {
    release();
    size     = size_;
    _data    = new Cell[size.prod()];
    pointers = new Cell*[size(1)];

    for (int j = 0; j < size(1); ++j) {
      pointers[j] = &_data[size(0) * j];
    }
  }

  void
  release() {
    if (_data) {
      delete _data;
    }

    if (pointers) {
      delete pointers;
    }
  }

  VectorDi size;
  Cell*    _data;
  Cell**   pointers;
};
template <typename Cell>
struct Pointers<Cell, 3> {
  typedef Eigen::Matrix<int, 3, 1> VectorDi;

  Pointers() : _data(0),
               pointers(0) {}

  Pointers(Pointers const& other) : size(other.size),
                                    _data(other._data),
                                    pointers(other.pointers) {}

  ~Pointers() {
    release();
  }

  Pointers&
  operator=(Pointers const& other) {
    size     = other.size;
    _data    = other._data;
    pointers = other.pointers;

    return *this;
  }

  Cell*
  getCell(VectorDi const& index) {
    return &pointers[index(2)][index(1)][index(0)];
  }

  void
  allocate(VectorDi const& size_) {
    release();
    size     = size_;
    _data    = new Cell[size.prod()];
    pointers = new Cell * *[size(2)];

    for (int k = 0; k < size(2); ++k) {
      pointers[k] = new Cell*[size(1)];

      for (int j = 0; j < size(1); ++j) {
        pointers[j][k] =
          &_data[size(1) * size(0) * k + size(0) * j];
      }
    }
  }

  void
  release() {
    if (_data) {
      delete _data;
    }

    if (pointers) {
      for (int k = 0; k < size(2); ++k) {
        delete pointers[k];
      }
      delete pointers;
    }
  }

  VectorDi size;
  Cell*    _data;
  Cell***  pointers;
};
}

template <typename TCell, int D>
class CellHeap {
public:
  typedef TCell                      Cell;
  typedef Eigen::Matrix<int, D, 1>   VectorDi;
  typedef Private::Pointers<Cell, D> Pointers;

public:
  CellHeap() {}

  CellHeap(CellHeap const& other) = delete;

  ~CellHeap() {}

  CellHeap&
  operator=(CellHeap const& other) = delete;

  void
  allocate(VectorDi const& size) {
    _pointers.allocate(size);
  }

  void
  release() {
    _pointers.release();
  }

  Cell*
  getCell(VectorDi const& index) {
    return _pointers.getCell(index);
  }

  Cell*
  operator()(VectorDi const& index) {
    return getCell(index);
  }

private:
  Pointers _pointers;
};
}
#endif
