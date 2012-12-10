#include "common_object.h"

void initial() 
{
  int x, y, pp;
  my_double rho,nu,gradP;
  my_double  usq, vsq, u2, v2;
  my_double ui, vi, sumsq, sumsq2, uv;
  my_double Ts;
  my_double beta, betas, gravity_y;
  pop p_eq;
  my_double fn,kn;

#ifdef TEMPERATURE_MELTING
  Ts = property.T_solid;
#endif

#ifdef FLUID
  nu = property.nu;
#ifdef FLUID_FORCING_POISEUILLE
  gradP=property.gradP;
#endif
#endif

#ifdef TEMPERATURE_BUOYANCY
  beta = property.beta_t;
  gravity_y = property.gravity_y;
#endif
#ifdef SALT
#ifdef SALT_BUOYANCY
  betas = property.beta_s;
  gravity_y = property.gravity_y;
#endif
#endif

  /* initialize random seed: */
  srand( time(NULL) );



    for (y=0; y<NY+2; y++)
      for (x=0; x<NX+2; x++){ 

	/* constant density */
  for (pp = 0; pp < 9; pp++)
	p[IDX(y,x)].p[pp] = wgt[pp];

#if defined(TEMPERATURE_BUOYANCY) && defined(SALT_BUOYANCY)
 for (pp = 0; pp < 9; pp++)
  p[IDX(y,x)].p[pp] = wgt[pp]* (exp(beta*gravity_y*(y-0.5)*( (property.T_bot-property.T_ref) - 0.5*(property.deltaT/NY)*(y-0.5) )/cs2 
                                  + betas*gravity_y*(y-0.5)*( (property.S_bot-property.S_ref) - 0.5*(property.deltaS/NY)*(y-0.5))/cs2 ));
#else
#ifdef TEMPERATURE_BUOYANCY
  for (pp = 0; pp < 9; pp++)
    /* hydrostatic pressure */
    // fprintf(stderr,"temperature : Starting with hydrostatic density profile.\n");
    p[IDX(y,x)].p[pp] = wgt[pp]* (exp(beta*gravity_y*(y-0.5)*( (property.T_bot-property.T_ref) - 0.5*(property.deltaT/NY)*(y-0.5) )/cs2 ));
    // p[IDX(y,x)].p[pp] = wgt[pp]*(1.0 +  beta*gravity_y*(y-0.5)*( (property.T_bot-property.T_ref) - 0.5*(property.deltaT/NY)*(y-0.5) )/cs2 );
#endif
#ifdef SALT_BUOYANCY
  for (pp = 0; pp < 9; pp++)
    /* hydrostatic pressure */
    // fprintf(stderr,"salt : Starting with hydrostatic density profile.\n");
    p[IDX(y,x)].p[pp] = wgt[pp]* (exp(betas*gravity_y*(y-0.5)*( (property.S_bot-property.S_ref) - 0.5*(property.deltaS/NY)*(y-0.5) )/cs2 ));
  //	p[IDX(y,x)].p[pp] += wgt[pp]*( betas*gravity_y*(y-0.5)*( (property.S_bot-property.S_ref) - 0.5*(property.deltaS/NY)*(y-0.5) )/cs2 );
#endif
#endif

#ifdef TEMPERATURE_INITIAL_CONSTANT
  for (pp = 0; pp < 9; pp++)
	p[IDX(y,x)].p[pp] = wgt[pp];
#endif

#ifdef FLUID
#ifdef FLUID_INITIAL_POISEUILLE
	/* Poiseuille profile v_x(y) = 1/2 * (gradP / nu) * y * (L_y - y) */
	/* v_x(max) = 1/8 * (gradP / nu) * L_y^2 */
	//p[y][x].p[pp] = cx[pp]*wgt[pp]*(0.5*gradP/nu) * ((my_double)y-0.5) * ((my_double)NY-((my_double)y-0.5));

	/* horizontal */	

     p[IDX(y,x)].p[1] = wgt[1]*(gradP/nu) * ((double)y-0.5) * ((double)NY-((double)y-0.5)); 
     p[IDX(y,x)].p[5] = wgt[5]*(gradP/nu) * ((double)y-0.5) * ((double)NY-((double)y-0.5));
     p[IDX(y,x)].p[8] = wgt[8]*(gradP/nu) * ((double)y-0.5) * ((double)NY-((double)y-0.5));
     p[IDX(y,x)].p[6] = -wgt[6]*(gradP/nu) * ((double)y-0.5) * ((double)NY-((double)y-0.5)); 
     p[IDX(y,x)].p[3] = -wgt[3]*(gradP/nu) * ((double)y-0.5) * ((double)NY-((double)y-0.5)); 
     p[IDX(y,x)].p[7] = -wgt[7]*(gradP/nu) * ((double)y-0.5) * ((double)NY-((double)y-0.5));
     
     p[IDX(y,x)].p[1] *= m(p[IDX(y,x)]);
     p[IDX(y,x)].p[5] *= m(p[IDX(y,x)]);
     p[IDX(y,x)].p[8] *= m(p[IDX(y,x)]);
     p[IDX(y,x)].p[3] *= m(p[IDX(y,x)]);
     p[IDX(y,x)].p[6] *= m(p[IDX(y,x)]);
     p[IDX(y,x)].p[7] *= m(p[IDX(y,x)]);
     	


#endif
#endif

#ifdef FLUID
#ifdef METHOD_FORCING_GUO
	force[IDX(y,x)].x = force[IDX(y,x)].y = 0.0;
#endif
      v[IDX(y,x)].vx = vx(p[IDX(y,x)]);
      v[IDX(y,x)].vy = vy(p[IDX(y,x)]);
      dens[IDX(y,x)] = m(p[IDX(y,x)]);
#else
      /* for better numerical accuracy */
      v[IDX(y,x)].vx = 0.0; 
      v[IDX(y,x)].vy = 0.0; 
      dens[IDX(y,x)] = 1.0; //density is assumed to be 1
#endif
      } /* end for on x, y*/   
  
#ifdef FLUID
#ifdef FLUID_INITIAL_PERTURBATION
/*********************************************                                                                                                    
 *         x  y                               *
 *  0    (+0,+0)         6    2    5          *
 *  1    (+1,+0)          \   |   /           *
 *  2    (+0,+1)           \  |  /            *
 *  3    (-1,+0)            \ | /             *
 *  4    ( 0,-1)      3 <---- 0 ----> 1       *
 *  5    (+1,+1)            / | \             *
 *  6    (-1,+1)           /  |  \            *
 *  7    (-1,-1)          /   |   \           *
 *  8    (+1,-1)         7    4    8          * 
 *                                            *
 *********************************************/
    fn=0.00000000001;
    kn=10.0;
    for (y=0; y<NY+2; y++)
      for (x=0; x<NX+2; x++){ 
	
	kn=sin(2.*3.14*y/NY)*5.0;
	fn=sin(2.*3.14*x/NX)*0.2;
	/* horizontal */	

	p[IDX(y,x)].p[1] += wgt[1]*fn*sin(kn*2.*3.14*y/NY); 
	p[IDX(y,x)].p[5] += wgt[5]*fn*sin(kn*2.*3.14*y/NY); 
	p[IDX(y,x)].p[8] += wgt[8]*fn*sin(kn*2.*3.14*y/NY); 
	p[IDX(y,x)].p[6] += -wgt[6]*fn*sin(kn*2.*3.14*y/NY); 
	p[IDX(y,x)].p[3] += -wgt[3]*fn*sin(kn*2.*3.14*y/NY); 
	p[IDX(y,x)].p[7] += -wgt[7]*fn*sin(kn*2.*3.14*y/NY); 
	
	/* veritcal */	

	p[IDX(y,x)].p[6] += wgt[6]*fn*sin(kn*2.*3.14*x/NX); 
	p[IDX(y,x)].p[2] += wgt[2]*fn*sin(kn*2.*3.14*x/NX);  
	p[IDX(y,x)].p[5] += wgt[5]*fn*sin(kn*2.*3.14*x/NX); 
	p[IDX(y,x)].p[7] += -wgt[7]*fn*sin(kn*2.*3.14*x/NX); 
	p[IDX(y,x)].p[4] += -wgt[4]*fn*sin(kn*2.*3.14*x/NX);  
	p[IDX(y,x)].p[8] += -wgt[8]*fn*sin(kn*2.*3.14*x/NX); 
	

      } /* end for on x, y*/  
#endif
#endif   

#ifdef TEMPERATURE
    for (y=0; y<NY+2; y++)
      for (x=0; x<NX+2; x++){ 
	/* linear temperature gradient */
	tt[IDX(y,x)] = ( (property.T_bot-property.T_ref) - (property.deltaT/NY)*(y-0.5) );
#ifdef TEMPERATURE_INITIAL_PERTURBATION
	if(x<NX/2){ tt[IDX(y,x)] += 1.e-4; }else{ tt[IDX(y,x)] -= 1.e-4; }
#endif
	/* on the populations */
	for (pp = 0; pp < 9; pp++) 
	  g[IDX(y,x)].p[pp] = wgt[pp]*tt[IDX(y,x)];

#ifdef TEMPERATURE_INITIAL_CONSTANT
        /* constant salinity */
        tt[IDX(y,x)] = property.T_top;
        /* on the populations */
        for (pp = 0; pp < 9; pp++)
          g[IDX(y,x)].p[pp] = wgt[pp]*tt[IDX(y,x)];
#endif
      }  
#endif

#ifdef SALT
    for (y=0; y<NY+2; y++)
      for (x=0; x<NX+2; x++){ 
#ifdef SALT_INITIAL_LINEAR
	/* linear salinity gradient */
	ss[IDX(y,x)] = ( (property.S_bot-property.S_ref) - (property.deltaS/NY)*(y-0.5) );
#ifdef SALT_INITIAL_PERTURBATION
	if(x<NX/2){ ss[IDX(y,x)] += 1.e-5;}else{ ss[IDX(y,x)] -= 1.e-5; } 
#endif
	/* on the populations */
	for (pp = 0; pp < 9; pp++)
	  s[IDX(y,x)].p[pp] = wgt[pp]*ss[IDX(y,x)];
#endif
#ifdef SALT_INITIAL_CONSTANT
        /* constant salinity */
        ss[IDX(y,x)] = property.S_top;
        /* on the populations */
        for (pp = 0; pp < 9; pp++)
          s[IDX(y,x)].p[pp] = wgt[pp]*ss[IDX(y,x)];
#endif

      }  
#endif

#ifdef TEMPERATURE_MELTING
  /* 1 is fluid , 0 is solid */
  for (y=0; y<NY+2; y++)
    for (x=0; x<NX+2; x++){ 
#ifdef TEMPERATURE_MELTING_INITIAL_SOLID
      ll[IDX(y,x)]=llold[IDX(y,x)]=0.0;
#ifdef TEMPERATURE_MELTING_INITIAL_SOLID_CAVITY
  for (y=0; y<NY+2; y++){
    for (x=0; x<NX+2; x++) 
      /* cube */
      // if(y<=20 && fabs((float)x-(float)NX/2)<=10)  ll[IDX(y,x)]=llold[IDX(y,x)]=1.0;
      /* half a circle */
      if( y*y + pow((float)x-(float)NX/2,2.0) <= 900)  ll[IDX(y,x)]=llold[IDX(y,x)]=1.0;
  }
#endif
#else
      ll[IDX(y,x)]=llold[IDX(y,x)]=1.0;
#endif
      //      if(y>NY/2) ll[IDX(y,x)]=llold[IDX(y,x)]=1.0;
    }  
#ifdef SALT
  /* 1 is fluid , 0 is solid */
  for (y=0; y<NY+2; y++){
    //fprintf(stdout,"%d %g %g %g\n",y, tt[IDX(y,x)], ss[IDX(y,x)], Ts);
    for (x=0; x<NX+2; x++){ 
      ll[IDX(y,x)]=llold[IDX(y,x)]=1.0;
      //Ts = property.liquidus_slope * ss[IDX(y,x)];
      // if( tt[IDX(y,x)] < Ts ){ ll[IDX(y,x)]=llold[IDX(y,x)]=0.0; }
    }  
  }
#endif 
#endif

#ifdef NON_NEWTON
  for (y=0; y<NY+2; y++)
    for (x=0; x<NX+2; x++){ 
      tau[IDX(y,x)]=property.tau;
    }   
#endif

#ifdef THERMIC_SOLID
  for (y=0; y<NY+2; y++)
    for (x=0; x<NX+2; x++){ 
      tauT[IDX(y,x)]=property.tau_t;
    }   
#endif

}
