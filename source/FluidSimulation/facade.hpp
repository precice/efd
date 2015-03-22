#pragma once

#include <memory>

namespace FsiSimulation {
namespace FluidSimulation {
class SimulationController;
class Configuration;

std::unique_ptr<SimulationController>
create_simulation_controller(Configuration* configuration);

}
}
