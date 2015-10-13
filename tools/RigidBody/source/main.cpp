#include "Application.hpp"

int
main(int argc, char** argv) {
  using Structure::Application;

  Application application(argc, argv);

  int status = application.parseArguments();

  if (status == 1) {
    return 0;
  }

  application.initialize();

  application.run();

  return 0;
}
