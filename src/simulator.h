#ifndef SIMULATOR_H
#define SIMULATOR_H
#include "input_handler.h"

class Simulator : public InputHandler
{
public:
  Simulator();
  virtual ~Simulator();
  virtual int Idle();
  virtual int PaintGL();
  class SimulatorInternal;
  SimulatorInternal* simulator_;
  DECLARE_SINGLETON_CLASS(Simulator)
};

#endif // SIMULATOR_H
