#include "Grid.hpp"

#include <Uni/Logging/macros>

namespace FsiSimulation {
namespace FluidSimulation {
template <typename T>
void
Grid<T>::
initialize(VectorDiType const& size, FactoryType const&  factory) {
  _factory = factory;

  this->BaseType::initialize(size,
                             VectorDiType::Zero(),
                             VectorDiType::Zero(),
                             _factory);

  VectorDiType indent(VectorDiType::Ones());

  for (int i = 0; i < Dimensions; ++i) {
    VectorDiType leftIndent(VectorDiType::Zero());
    leftIndent(i) = 0;
    VectorDiType rightIndent(VectorDiType::Zero());
    rightIndent(i) = size(i) - 1;

    boundaries[i][0].initialize(size,
                                leftIndent,
                                rightIndent,
                                _factory);

    boundaries[i][1].initialize(size,
                                rightIndent,
                                leftIndent,
                                _factory);
    leftIndent     = indent;
    leftIndent(i)  = 0;
    rightIndent    = indent;
    rightIndent(i) = size(i) - rightIndent(i);

    indentedBoundaries[i][0].initialize(size,
                                        leftIndent,
                                        rightIndent,
                                        _factory);

    indentedBoundaries[i][1].initialize(size,
                                        rightIndent,
                                        leftIndent,
                                        _factory);
  }

  innerGrid.initialize(size, indent, indent, _factory);
}

template <typename T>
std::string
Grid<T>::
toString() const {
  INFO << "Grid size: "
       << this->size().transpose() << "\n"
       << "Grid left indent: "
       << this->leftIndent().transpose() << "\n"
       << "Grid right indent: "
       << this->rightIndent().transpose() << "\n"
       << "Inner grid size: "
       << this->innerGrid.size().transpose() << "\n"
       << "Inner grid left indent: "
       << this->innerGrid.leftIndent().transpose() << "\n"
       << "Inner grid right indent: "
       << this->innerGrid.rightIndent().transpose() << "\n";

  for (unsigned d = 0; d < Dimensions; ++d) {
    INFO << 2 * d << " grid size: "
         << this->boundaries[d][0].size().transpose() << "\n"
         << 2 * d << " grid left indent: "
         << this->boundaries[d][0].leftIndent().transpose() << "\n"
         << 2 * d << " grid right indent: "
         << this->boundaries[d][0].rightIndent().transpose() << "\n"
         << 2 * d + 1 << " grid size: "
         << this->boundaries[d][1].size().transpose() << "\n"
         << 2 * d + 1 << " grid left indent: "
         << this->boundaries[d][1].leftIndent().transpose() << "\n"
         << 2 * d + 1 << " grid right indent: "
         << this->boundaries[d][1].rightIndent().transpose() << "\n";
  }

  return "";
}

Fluid_InstantiateExternTemplates(Grid);
}
}
