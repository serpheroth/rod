#include <memory>
#include "simulator.h"
#include "pbd_elastic_rod.h"

class Simulator::SimulatorInternal
{
public:
  SimulatorInternal()
    : rod_(new PBDElasticRod(rod::rod_v_num, rod::rod_length))
  {
  }

  ~SimulatorInternal()
  {
    delete rod_;
  }
  PBDElasticRod* rod_;
};

Simulator::Simulator()
  : InputHandler(0)
  , simulator_(new SimulatorInternal)
{
}

int Simulator::Idle()
{
  for (int i = 0; i < global::simulation_step_per_idle; ++i) {
    simulator_->rod_->Simulate(global::time_step);
  }
  return -1;
}

int Simulator::PaintGL()
{
  simulator_->rod_->Render();
  return -1;
}

Simulator::~Simulator()
{
  delete simulator_;
}
