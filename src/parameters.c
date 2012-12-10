#include "common_object.h"

void design_lb(){

  /* weights */
  my_double rt0=(4.0 /  9.0);
  my_double rt1=(1.0 /  9.0);
  my_double rt2=(1.0 / 36.0);

  wgt[0] = rt0; 
  wgt[1] = rt1;
  wgt[2] = rt1;
  wgt[3] = rt1;
  wgt[4] = rt1;
  wgt[5] = rt2;
  wgt[6] = rt2;
  wgt[7] = rt2;
  wgt[8] = rt2;

  /* lattice velocities */
  cx[0] =  0.0;
  cx[1] =  1.0;  
  cx[2] =  0.0;
  cx[3] = -1.0;
  cx[4] =  0.0;  
  cx[5] =  1.0;
  cx[6] = -1.0;
  cx[7] = -1.0;  
  cx[8] =  1.0;

  cy[0] =  0.0;
  cy[1] =  0.0;  
  cy[2] =  1.0;
  cy[3] =  0.0;
  cy[4] = -1.0;  
  cy[5] =  1.0;
  cy[6] =  1.0;
  cy[7] = -1.0;  
  cy[8] = -1.0;

  /* directions on the lattice */
  /*
  dirp = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  invp = {0, 3, 4, 1, 2, 7, 8, 5, 6};
  */
  /* speed of sound constants */
extern my_double cs, cs2 , cs4 , cs22 , cssq;
extern my_double invcs, invcs2, invcs4;

  cs=1.0/sqrt(3.0);   invcs = 1.0/cs;
  cs2=(1.0 / 3.0);    invcs2 = 1.0/cs2;
  cs4=(1.0 / 9.0);    invcs4 = 1.0/cs4;
  twocs2=2.0*cs2; invtwocs2 = 1.0/twocs2;
  twocs4=2.0*cs4; invtwocs4 = 1.0/twocs4;
}


my_double read_parameter(char * variable){

  char fnamein[256], fnameout[256];
  char name[256] = "NULL";
  double val;
  FILE *fin, *fout;
  int i;
  int cmp=1;

  sprintf(fnamein,"param.in");
  sprintf(fnameout,"param.out");
  fin = fopen(fnamein,"r");
  if(fin == NULL){
	 fprintf(stderr,"Error message -> %s file is missing! Exit.\n",fnamein);
	 fclose(fin);
	 exit(0);
  }

  while(cmp!=0){  
    val=0.0;
    fscanf(fin,"%s %lf",&name,&val);
    cmp=strcmp(name,variable);
       if(cmp==0){
	 fout = fopen(fnameout,"a");
	 fprintf(fout,"%s %g\n",name, val);
	 fclose(fout);
	 return (my_double)val;
       }
       if(feof(fin)){ 
	 fprintf(stderr,"Error message -> %s not found. End of File. Exit.\n",variable);
	 fclose(fin);
	 exit(0);
       }    
  }
  fclose(fin);
}

void assign_parameters(){
  char name[256];

  sprintf(OutDir,"RUN");
  mkdir(OutDir, S_IWUSR|S_IXUSR|S_IRUSR);
  fprintf(stderr,"OutDir is %s\n",OutDir);   

  remove("param.out");
  /* read parameters from file */
  sprintf(name,"NX");
  property.NX = NX = (int)read_parameter(name);
  sprintf(name,"NY");
  property.NY = NY = (int)read_parameter(name);
  fprintf(stderr,"System Size:\nNX %d \nNY %d\n", NX , NY);

  sprintf(name,"max_step");
  max_step = (int)read_parameter(name); 
  fprintf(stderr,"Total time steps: %d\n",max_step);
  
  sprintf(name,"time_dump_field");
  time_dump_field = (int)read_parameter(name);
  fprintf(stderr,"Time dump fields: %d\n",time_dump_field);

#ifdef FLUID
  fprintf(stderr,"YES <- FLUID\n");
  sprintf(name,"tau");
  property.tau = read_parameter(name);
  fprintf(stderr,"Properties:\ntau %g\n",(double)property.tau);
  property.nu = (property.tau - 0.5)/3.0;
  fprintf(stderr,"viscosity %g\n",(double)property.nu);
#ifdef FLUID_FORCING_POISEUILLE
  fprintf(stderr,"YES <- FLUID_FORCING_POISEUILLE\n");
  sprintf(name,"gradP");
  property.gradP = read_parameter(name);
  fprintf(stderr,"Properties:\ngradP %g\n",(double)property.gradP);
  dimensionless.Reynolds = ((2./3.)*(double)property.gradP/(2.*(double)property.nu)*(double)property.NY*(double)property.NY/4.0)*(double)property.NY/(double)property.nu;
  /* NOTE :density shall be corrected */
  fprintf(stderr,"Reynolds %g\n",(double)dimensionless.Reynolds);
#endif
#endif

#ifdef TEMPERATURE
  fprintf(stderr,"YES <- TEMPERATURE\n");
  sprintf(name,"tau_t");
  property.tau_t = read_parameter(name);
  fprintf(stderr,"Properties:\ntau_t %g\n",(double)property.tau_t);
  property.kappa_t = (property.tau_t - 0.5)/3.0;
  fprintf(stderr,"thermal diffusivity %g\n",(double)property.kappa_t);
  sprintf(name,"T_bot");
  property.T_bot = read_parameter(name);
  sprintf(name,"T_top");
  property.T_top = read_parameter(name);
  sprintf(name,"T_ref");
  property.T_ref = read_parameter(name);
  property.deltaT = property.T_bot-property.T_top;
  fprintf(stderr,"T_bot %g , T_top %g , deltaT %g\n",(double)property.T_bot, (double)property.T_top, (double)property.deltaT);
#ifdef TEMPERATURE_BUOYANCY
  fprintf(stderr,"YES <- TEMPERATURE_BUOYANCY\n");
  sprintf(name,"beta_t");
  property.beta_t = read_parameter(name);
  fprintf(stderr,"linear volume expansion coefficient %g\n",(double)property.beta_t);
  sprintf(name,"beta2_t");
  property.beta2_t = read_parameter(name);
  fprintf(stderr,"quadratic volume expansion coefficient %g\n",(double)property.beta2_t);
  sprintf(name,"gravity_y");
  property.gravity_y = read_parameter(name);
  fprintf(stderr,"gravity_y %g\n",(double)property.gravity_y);
#endif 
#ifdef TEMPERATURE_MELTING
  fprintf(stderr,"YES <- TEMPERATURE_MELTING\n");
  sprintf(name,"T_solid");
  property.T_solid = read_parameter(name);
  sprintf(name,"specific_heat");
  property.specific_heat = read_parameter(name);
  sprintf(name,"latent_heat");
  property.latent_heat = read_parameter(name);
  sprintf(name,"liquidus_slope");
  property.liquidus_slope = read_parameter(name);
  fprintf(stderr,"Temperature solid %g\n",(double)property.T_solid);
  fprintf(stderr,"Thermal capacity at constant pressure %g\n",(double)property.specific_heat);
  fprintf(stderr,"Latent_heat %g\n",(double)property.latent_heat);
  fprintf(stderr,"Liquidus_slope %g\n",(double)property.liquidus_slope);
#endif 
#endif

#ifdef SALT
  fprintf(stderr,"YES <- SALT\n");
  sprintf(name,"tau_s");
  property.tau_s = read_parameter(name);
  fprintf(stderr,"Properties:\ntau_s %g\n",(double)property.tau_s);
  property.kappa_s = (property.tau_s - 0.5)/3.0;
  fprintf(stderr,"salt diffusivity %g\n",(double)property.kappa_s);
  sprintf(name,"S_bot");
  property.S_bot = read_parameter(name);
  sprintf(name,"S_top");
  property.S_top = read_parameter(name);
  sprintf(name,"S_ref");
  property.S_ref = read_parameter(name);
  property.deltaS = property.S_bot-property.S_top;
  fprintf(stderr,"S_bot %g , S_top %g , deltaS %g\n",(double)property.S_bot, (double)property.S_top, (double)property.deltaS);
#ifdef SALT_BUOYANCY
  fprintf(stderr,"YES <- SALT_BUOYANCY\n");
  sprintf(name,"beta_s");
  property.beta_s = read_parameter(name);
  fprintf(stderr,"linear volume expansion coefficient %g\n",(double)property.beta_t);
#ifndef TEMPERATURE_BUOYANCY
  sprintf(name,"gravity_y");
  property.gravity_y = read_parameter(name);
  fprintf(stderr,"gravity_y %g\n",(double)property.gravity_y);
#endif
#endif 
#endif

#ifdef TEMPERATURE 
#ifdef FLUID
 fprintf(stderr,"\n -- DIMENSIONLESS CONTROL PARAMETERS -- \n");
  dimensionless.Prandtl = property.nu/property.kappa_t;
  fprintf(stderr,"Prandtl %g\n",(double)dimensionless.Prandtl);
#ifdef TEMPERATURE_BUOYANCY
  dimensionless.Rayleigh_t = property.beta_t*property.gravity_y*property.deltaT*pow((double)property.NY,3.0)/(property.nu*property.kappa_t);
  fprintf(stderr,"Thermal Rayleigh %g\n",(double)dimensionless.Rayleigh_t);
#endif
#endif
#ifdef TEMPERATURE_MELTING
  dimensionless.Stefan= property.deltaT*property.specific_heat/property.latent_heat;
  fprintf(stderr,"Stefan %g\n",(double)dimensionless.Stefan);
#endif
#endif
#ifdef SALT
#ifdef FLUID
  dimensionless.Schmidt = property.nu/property.kappa_s;
  fprintf(stderr,"Schmidt %g\n",(double)dimensionless.Schmidt);
#endif
#endif
#ifdef SALT && TEMPERATURE
  dimensionless.Lewis = property.kappa_t/property.kappa_s;
  fprintf(stderr,"Lewis %g\n",(double)dimensionless.Lewis);
#endif
#ifdef SALT
#ifdef FLUID
#ifdef SALT_BUOYANCY
  dimensionless.Rayleigh_s = property.beta_s*property.gravity_y*property.deltaS*pow((double)property.NY,3.0)/(property.nu*property.kappa_s);
  fprintf(stderr,"Solutal Rayleigh %g\n",(double)dimensionless.Rayleigh_s);
#endif
#endif
#endif
 fprintf(stderr,"\n");

 fprintf(stderr,"Size of float %d\n",sizeof(float));
 fprintf(stderr,"Size of double %d\n",sizeof(double));
 fprintf(stderr,"Size of long double %d\n",sizeof(long double));
 fprintf(stderr,"Size of my_double %d\n",sizeof(my_double));
 fprintf(stderr,"\n");
}

void allocate_fields(){
  //#ifdef FLUID
  p  = (pop*) malloc(sizeof(pop)*(NX+2)*(NY+2)); 
 if(p == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}

#ifdef METHOD_STEPPING_AB2
  p_old  = (pop*) malloc(sizeof(pop)*(NX+2)*(NY+2)); 
 if(p_old == NULL){ fprintf(stderr,"Not enough memory to allocate p_old\n"); exit(-1);}
#endif

  buffer  = (pop*) malloc(sizeof(pop)*(NX+2)*2); 
 if(buffer == NULL){ fprintf(stderr,"Not enough memory to allocate buffer\n"); exit(-1);}

  v  = (velocity*) malloc(sizeof(velocity)*(NX+2)*(NY+2)); 
 if(v == NULL){ fprintf(stderr,"Not enough memory to allocate v\n"); exit(-1);}

  vold  = (velocity*) malloc(sizeof(velocity)*(NX+2)*(NY+2)); 
 if(vold == NULL){ fprintf(stderr,"Not enough memory to allocate vold\n"); exit(-1);}

  dens  = (my_double*) malloc(sizeof(my_double)*(NX+2)*(NY+2)); 
 if(dens == NULL){ fprintf(stderr,"Not enough memory to allocate dens\n"); exit(-1);}

  force  = (vector*) malloc(sizeof(vector)*(NX+2)*(NY+2)); 
 if(force == NULL){ fprintf(stderr,"Not enough memory to allocate force\n"); exit(-1);}

#ifdef FLUID_RHEOLOGY
  tau  = (my_double*) malloc(sizeof(my_double)*(NX+2)*(NY+2)); 
 if(tau == NULL){ fprintf(stderr,"Not enough memory to allocate tau\n"); exit(-1);}
#endif
 //#endif 
#ifdef FLUID_GRADIENT
  gradv  = (tensor*) malloc(sizeof(vector)*(NX+2)*(NY+2)); 
 if(gradv == NULL){ fprintf(stderr,"Not enough memory to allocate gradv\n"); exit(-1);}
#endif

#ifdef FLUID_RHEOLOGY_POLYMER
  conf  = (tensor*) malloc(sizeof(vector)*(NX+2)*(NY+2)); 
 if(conf == NULL){ fprintf(stderr,"Not enough memory to allocate conf\n"); exit(-1);}
#endif
#ifdef FLUID_POROSITY
 porosity  = (my_double*) malloc(sizeof(my_double)*(NX+2)*(NY+2));
 if(porosity == NULL){ fprintf(stderr,"Not enough memory to allocate porosity\n"); exit(-1);}
#endif

#ifdef TEMPERATURE
  g  = (pop*) malloc(sizeof(pop)*(NX+2)*(NY+2)); 
 if(g == NULL){ fprintf(stderr,"Not enough memory to allocate g\n"); exit(-1);}

#ifdef METHOD_STEPPING_AB2
  g_old  = (pop*) malloc(sizeof(pop)*(NX+2)*(NY+2)); 
 if(g_old == NULL){ fprintf(stderr,"Not enough memory to allocate g_old\n"); exit(-1);}
#endif

  tt  = (my_double*) malloc(sizeof(my_double)*(NX+2)*(NY+2)); 
 if(tt == NULL){ fprintf(stderr,"Not enough memory to allocate tt\n"); exit(-1);}

  ttold  = (my_double*) malloc(sizeof(my_double)*(NX+2)*(NY+2)); 
 if(ttold == NULL){ fprintf(stderr,"Not enough memory to allocate ttold\n"); exit(-1);}

#ifdef TEMPERATURE_MELTING
  ll  = (my_double*) malloc(sizeof(my_double)*(NX+2)*(NY+2)); 
 if(ll == NULL){ fprintf(stderr,"Not enough memory to allocate ll\n"); exit(-1);}

  llold  = (my_double*) malloc(sizeof(my_double)*(NX+2)*(NY+2)); 
 if(llold == NULL){ fprintf(stderr,"Not enough memory to allocate llold\n"); exit(-1);}

  hh  = (my_double*) malloc(sizeof(my_double)*(NX+2)*(NY+2)); 
 if(hh == NULL){ fprintf(stderr,"Not enough memory to allocate hh\n"); exit(-1);}
#endif
#endif

#ifdef SALT
  s  = (pop*) malloc(sizeof(pop)*(NX+2)*(NY+2)); 
 if(s == NULL){ fprintf(stderr,"Not enough memory to allocate g\n"); exit(-1);}

#ifdef METHOD_STEPPING_AB2
  s_old  = (pop*) malloc(sizeof(pop)*(NX+2)*(NY+2)); 
 if(s_old == NULL){ fprintf(stderr,"Not enough memory to allocate s_old\n"); exit(-1);}
#endif

  ss  = (my_double*) malloc(sizeof(my_double)*(NX+2)*(NY+2)); 
 if(ss == NULL){ fprintf(stderr,"Not enough memory to allocate ss\n"); exit(-1);}

  ssold  = (my_double*) malloc(sizeof(my_double)*(NX+2)*(NY+2)); 
 if(ssold == NULL){ fprintf(stderr,"Not enough memory to allocate ssold\n"); exit(-1);}
#endif
}


//int IDX(int j, int i){
//return (j*(NX+2)+i);
  //return (i*(NY+2)+j);
//}




