#include "define.h"

void design_lb();
double read_parameter();
void assign_parameters();
void allocate_fields();
void initial();
//int IDX();

#ifdef FLUID
pop equilibrium();
void displace(); 
void collide();
void bc();
void hydro_fields();
void check_mass();
#ifdef FLUID_FORCING_POISEUILLE
void poiseuille_forc();
#endif
#endif

#ifdef METHOD_FORCING_GUO
void apply_forcing();
#endif

#ifdef TEMPERATURE
void bct();
void nusselt();
void collidet();
void displacet();
#ifdef TEMPERATURE_BUOYANCY
void buoyancy();
#endif
#endif
