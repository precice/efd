#ifndef FsiSimulation_MySimulation_hpp
#define FsiSimulation_MySimulation_hpp

namespace FsiSimulation {
class MySimulation {
public:
  virtual
  ~MySimulation() {}

  virtual void
  initialize() = 0;

  virtual void
  iterate() = 0;
};
}
#endif
