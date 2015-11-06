#include "Application.hpp"

int
main(int argc, char* argv[]) {
  using namespace Fluid::EntryPoint;

  Application application(argc, argv);

  auto result = application.parseArguments();

  if (result == 1) {
    return 0;
  }

  if (result != 0) {
    return result;
  }

  application.initialize();

  application.run();

  application.release();

  return 0;
}
