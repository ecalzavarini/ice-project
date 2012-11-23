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

#if defined(FLUID) || defined(TEMPERATURE)
     nusselt(itime,0);  /* compute heat flow and other global averages : 1st part */
     //nusselt(itime,1);
#endif

#ifdef NON_NEWTON
    relaxation_time(); /* compute local relaxation time */
#endif

    /* Collisions */
#ifdef FLUID
#ifdef METHOD_STEPPING_AB2
    collide(p,p_old,0);
#else
   collide(p,0);
#endif
#endif

#ifdef TEMPERATURE
#ifdef METHOD_STEPPING_AB2
    collide(g,g_old,1);
#else
   collide(g,1);
#endif
#endif

#ifdef SALT
#ifdef METHOD_STEPPING_AB2
    collide(s,s_old,2);
#else
   collide(s,2);
#endif
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

#ifdef FLUID_RHEOLOGY_POLYMER
    polymer_conformation_equation();
#endif

#ifdef FLUID
#ifdef METHOD_FORCING_GUO
    apply_forcing();
#endif
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

#if defined(FLUID) || defined(TEMPERATURE)
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
