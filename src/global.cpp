#include "global.h"
#include "macro_constant.h"
#include "config_file.h"


ConfigFile conf(DATA_DIRECTORY "global.conf");
namespace global
{
OpenGLQt* gl = NULL;
const char kPWD[] = DATA_DIRECTORY;
const char kDataDir[] = DATA_DIRECTORY;
std::string data_directory(DATA_DIRECTORY);

int simulation_step_per_idle = conf.Get<int>("simulation_step_per_idle");
int simulate = conf.Get<int>("simulate");
int pause_per_frame = conf.Get<int>("pause_per_frame");
real time_step = conf.Get<real>("time_step");
real* gravity = conf.Get<real*>("gravity");

}

namespace rod
{
int pbd_iteration = conf.Get<int>("pbd_iteration");
real* bending_twising_stiffness = conf.Get<real*>("bending_twisting_stiffness");
int rod_v_num = conf.Get<int>("rod_v_num");
real rod_length = conf.Get<real>("rod_length");
}


