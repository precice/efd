#pragma once

#include <memory>

namespace Fluid {
namespace Simulation {
class SimulationController;
class Configuration;

std::unique_ptr<SimulationController>
create_simulation_controller(Configuration* configuration);

}
}
