#include "typedef.h"
#include "common_main.h"

int main(int argc, char** argv){ 

  design_lb();  
  assign_parameters();
  allocate_fields();

  /* initialize the fields */
  initial();  

  print_fields(-100);
  /* enforce boundary conditions */

#ifdef FLUID
  bc();
#endif

#ifdef TEMPERATURE 
  bct();
#endif

#ifdef SALT
  bcs();
#endif 


  for (itime=0; itime<max_step; itime++) {
    if(itime%100==0) fprintf(stderr,"time step %d\n",itime);
  
    /* compute hydrodynamic fields an print them */
    hydro_fields(itime);      
    print_fields(itime);    
    check_mass(itime);

#ifdef TEMPERATURE
     nusselt(itime,0);  /* compute heat flow and other global averages : 1st part */
#endif

#ifdef NON_NEWTON
    relaxation_time(); /* compute local relaxation time */
#endif

    /* Collisions */
#ifdef FLUID
    collide(p,0);
#endif

#ifdef TEMPERATURE
    collide(g,1);
#endif

#ifdef SALT
    collide(s,2);
#endif

    /* Forcings */
#ifdef TEMPERATURE_BUOYANCY
     buoyancy();
#endif 

#ifdef SALT_BUOYANCY
    buoyancys();
#endif

#ifdef FLUID_FORCING_POISEUILLE
    poiseuille_forc();
#endif

#ifdef TEMPERATURE_MELTING
    melting(itime);
#endif

#ifdef METHOD_FORCING_GUO
    apply_forcing();
#endif


    /*  Boundary conditions */  
#ifdef FLUID   
      bc();
#endif

#ifdef TEMPERATURE
      bct();
#endif

#ifdef SALT
      bcs();
#endif

#ifdef TEMPERATURE
      nusselt(itime,1); 
#endif

    /*  Streaming */
#ifdef FLUID
    displace(p);
#endif

#ifdef TEMPERATURE
    displace(g);
#endif

#ifdef SALT
    displace(s);
#endif

  } /* end loop for on itime */

  return(0);
}
