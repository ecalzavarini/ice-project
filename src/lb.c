#include "common_object.h"

void displace(pop* p) 
{
  int x, y, idx;
  int xm, xp, ym, yp;
  int plane = 0;

  for (y=1; y<NY+1; y++) {
    ym = y-1;
    yp = y+1;
    for (x=1; x<NX+1; x++) {
      xm = x-1;
      xp = x+1;
      idx=IDXbuffer(plane,x);

      buffer[idx].p[0] = p[IDX(y,x)  ].p[0]; 
      buffer[idx].p[1] = p[IDX(y,xm) ].p[1]; 
      buffer[idx].p[5] = p[IDX(ym,xm)].p[5]; 
      buffer[idx].p[2] = p[IDX(ym,x) ].p[2]; 
      buffer[idx].p[6] = p[IDX(ym,xp)].p[6]; 
      buffer[idx].p[3] = p[IDX(y,xp) ].p[3]; 
      buffer[idx].p[7] = p[IDX(yp,xp)].p[7]; 
      buffer[idx].p[4] = p[IDX(yp,x) ].p[4]; 
      buffer[idx].p[8] = p[IDX(yp,xm)].p[8]; 
    }  
    
    if (plane) plane=0; else plane=1;
    for (x=1; x<NX+1; x++) 
      p[IDX(y-1,x)] = buffer[IDXbuffer(plane,x)];
  }

  if (plane) plane=0; else plane=1;
  for (x=1; x<NX+1; x++) 
    p[IDX(NY,x)] = buffer[idx=IDXbuffer(plane,x)];
}



pop equilibrium(pop *f, int y, int x) {
  int pp;
  my_double ux, uy;
  my_double rho,rhof;
  my_double cu, u2;
  pop f_eq;

      rhof = m(f[IDX(y,x)]);
      //rho = m(p[IDX(y,x)]);

      ux = v[IDX(y,x)].vx;
      uy = v[IDX(y,x)].vy;

      u2 = (ux*ux +  uy*uy);
  
      /* equilibrium distribution */
      for (pp=0; pp<9; pp++){
	cu = (cx[pp]*ux + cy[pp]*uy);
	//	f_eq.p[pp] = rhof * wgt[pp] * (1.0 + 3.0*cu  + 4.5*cu*cu - 1.5*u2 );
	f_eq.p[pp] = rhof * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*cu*cu - invtwocs2*u2 );
      }

      return f_eq;    
}

#ifdef FLUID_POROSITY
pop equilibrium_porosity(pop *f, int y, int x) {
  int pp;
  my_double ux, uy;
  my_double rho,rhof;
  my_double cu, u2 , inveps;
  pop f_eq;

  rhof = m(f[IDX(y,x)]);
  //rho = m(p[IDX(y,x)]);                                                                                                 

  ux = v[IDX(y,x)].vx;
  uy = v[IDX(y,x)].vy;

  u2 = (ux*ux +  uy*uy);

  inveps= 1./porosity[IDX(y,x)];

  /* equilibrium distribution */
  for (pp=0; pp<9; pp++){
    cu = (cx[pp]*ux + cy[pp]*uy);
    //      f_eq.p[pp] = rhof * wgt[pp] * (1.0 + 3.0*cu  + 4.5*cu*cu - 1.5*u2 );                                          
    f_eq.p[pp] = rhof * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*inveps*cu*cu - invtwocs2*inveps*u2 );
  }

  return f_eq;
}
#endif



pop equilibrium_opposite_velocity(pop *f, int y, int x) {
  int pp;
  my_double ux, uy;
  my_double rho,rhof;
  my_double cu, u2;
  pop f_eq;

      rhof = m(f[IDX(y,x)]);
      //rho = m(p[IDX(y,x)]);

      ux = v[IDX(y,x)].vx;
      uy = v[IDX(y,x)].vy;

      u2 = (ux*ux +  uy*uy);
  
      /* equilibrium distribution */
      for (pp=0; pp<9; pp++){
	cu = -(cx[pp]*ux + cy[pp]*uy);
	//	f_eq.p[pp] = rhof * wgt[pp] * (1.0 + 3.0*cu  + 4.5*cu*cu - 1.5*u2 );
	f_eq.p[pp] = rhof * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*cu*cu - invtwocs2*u2 );
      }

      return f_eq;    
}


#ifdef METHOD_STEPPING_AB2
void collide(pop* f,pop* f_old,int n){
#else
void collide(pop* f,int n){
#endif

  int x, y, pp;
  my_double invtau , one_minus_invtau;
  pop f_eq;

switch (n) {
#ifdef FLUID
 case 0:  /* This is fluid */
   invtau = 1.0 / property.tau;   
  break;
#endif
#ifdef TEMPERATURE
 case 1:  /* This is temperature */
   invtau = 1.0 / property.tau_t;  
  break;
#endif
#ifdef SALT
 case 2:  /* This is salt */
   invtau = 1.0 / property.tau_s;  
  break;
#endif
 } 
  
 one_minus_invtau = (1.0 - invtau);

  for (y=1; y<NY+1; y++) {
    for (x=1; x<NX+1; x++) {
      
      f_eq=equilibrium(f,y,x);
      /* collision */
      for (pp=0; pp<9; pp++)
#ifdef METHOD_STEPPING_AB2
	/* original */	
    //f[IDX(y,x)].p[pp] += -1.5*invtau * (f[IDX(y,x)].p[pp] - f_eq.p[pp]) + 0.5*invtau * (f_old[IDX(y,x)].p[pp] - f_eq.p[pp]);
	/* compact */
      f[IDX(y,x)].p[pp] += -0.5*invtau*(3.0*f[IDX(y,x)].p[pp]  - f_old[IDX(y,x)].p[pp] + 2.0*f_eq.p[pp]) ;
      /* store old */
      f_old[IDX(y,x)].p[pp] = f[IDX(y,x)].p[pp];
#else
      /* Euler Stepping */
      /* original */
	//	f[IDX(y,x)].p[pp] = f[IDX(y,x)].p[pp] - invtau * (f[IDX(y,x)].p[pp] - f_eq.p[pp]);
      /* compact */
	f[IDX(y,x)].p[pp] =  one_minus_invtau * f[IDX(y,x)].p[pp] + invtau * f_eq.p[pp];
#endif
    }  
  }

}

#ifdef FLUID
#ifdef METHOD_FORCING_GUO
/* Discrete lattice effects on the forcing term in the lattice Boltzmann method */
/* Zhaoli Guo, Chuguang Zheng, and Baochang Shi */
/* PHYSICAL REVIEW E, VOLUME 65, 046308 */
void apply_forcing() 
{
  my_double coeff,invtau;
  int x, y, pp, idx;
  my_double ux , uy , cu , dx, dy ;
  pop popForce;

   invtau = 1.0 / property.tau;  
    coeff = (1.0 - 0.5*invtau);

  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){
      idx=IDX(y,x);

        ux = v[idx].vx;
        uy = v[idx].vy;

      for (pp=0; pp<9; pp++){
	cu = (cx[pp]*ux + cy[pp]*uy);
	dx = (cx[pp]-ux)*invcs2 + cx[pp]*cu*invcs4;
	dy = (cy[pp]-uy)*invcs2 + cy[pp]*cu*invcs4;

	popForce.p[pp] = coeff * dens[idx] *wgt[pp]*(dx*force[idx].x  + dy*force[idx].y);
	p[idx].p[pp] += popForce.p[pp];
      }

    }
}
#endif
#endif



/**************************************************/
void hydro_fields(int i){
  int x,y;
  my_double tmpx, tmpy;
  FILE *ferr;
  my_double error;

  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++) {

#ifdef FLUID
      dens[IDX(y,x)] = m(p[IDX(y,x)]);
      #ifdef METHOD_FORCING_GUO
      v[IDX(y,x)].vx = ( vx(p[IDX(y,x)]) + 0.5*force[IDX(y,x)].x )/dens[IDX(y,x)];
      v[IDX(y,x)].vy = ( vy(p[IDX(y,x)]) + 0.5*force[IDX(y,x)].y )/dens[IDX(y,x)];
     /* set to zero after computing velocity */
      force[IDX(y,x)].x = force[IDX(y,x)].y = 0.0;
      #else
      v[IDX(y,x)].vx = vx(p[IDX(y,x)])/dens[IDX(y,x)];
      v[IDX(y,x)].vy = vy(p[IDX(y,x)])/dens[IDX(y,x)];
      #endif
#endif

#ifdef FLUID_POROSITY
      /* following Guo & Zhao , PHYSICAL REVIEW E 66, 036304 (2002) */
#ifdef METHOD_FORCING_GUO
      eps=porosity[IDX(y,x)];
      vx_temp = ( vx(p[IDX(y,x)]) + 0.5*eps*force[IDX(y,x)].x )/dens[IDX(y,x)];
      vy_temp = ( vy(p[IDX(y,x)]) + 0.5*eps*force[IDX(y,x)].y )/dens[IDX(y,x)];
      v_amp = sqrt(vx_temp*vx_temp + vy_temp*vy_temp); 
      c0 = 0.5 + 0.25*porperty.nu*(1-eps*eps)/(eps*eps);
      c1 = 0.0; 	
      v[IDX(y,x)].vx = vx_temp/(c0+sqrt(c0*c0 + c1*v_amp));
      v[IDX(y,x)].vy = vy_temp/(c0+sqrt(c0*c0 + c1*v_amp));
      /* set to zero after computing velocity */
      force[IDX(y,x)].x = force[IDX(y,x)].y = 0.0;
#endif
#endif

#ifdef TEMPERATURE
      tt[IDX(y,x)]    =  t(g[IDX(y,x)]);
#endif

#ifdef SALT
      ss[IDX(y,x)]    =  t(s[IDX(y,x)]);
#endif
    }

  /* check thermalization */
  if (i%500 == 0) {
    ferr  = fopen("error.dat","a");
    error = 0.0;
    for (y=1; y<NY+1; y++) 
      for (x=1; x<NX+1; x++) {
	tmpx = v[IDX(y,x)].vx - vold[IDX(y,x)].vx;
	tmpy = v[IDX(y,x)].vy - vold[IDX(y,x)].vy;  
	error += (tmpx * tmpx + tmpy * tmpy);
      }
    fprintf(ferr,"%d %g\n",i,error);
    fflush(ferr);

    if ((error < 10e-11) && (i!=0)) {
      fprintf(stderr,"Run termalized\n");
      fclose(ferr);
    }
  }/* if */

  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++) {
      vold[IDX(y,x)].vx = v[IDX(y,x)].vx;
      vold[IDX(y,x)].vy = v[IDX(y,x)].vy;
#ifdef TEMPERATURE
      ttold[IDX(y,x)] = tt[IDX(y,x)];
#endif
    }
}


/*************************************/
void check_mass(int i){
  int x, y;
  my_double mass1, mass2, mass3; 
  FILE *fmass;

  fmass = fopen("mass.dat","a");
  mass1 = mass2 = mass3 = 0.0;

  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++){ 
      mass1 += dens[IDX(y,x)];
    }
 fprintf(fmass,"%d %g ",i,mass1);

#ifdef TEMPERATURE
  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++){ 
      mass2 += tt[IDX(y,x)];
    }
 fprintf(fmass,"%g ",mass2); 
#endif

#ifdef SALT
  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++){ 
      mass3 += ss[IDX(y,x)];
    }
 fprintf(fmass,"%g ",mass3); 
#endif
  
 fprintf(fmass,"\n"); 

  fflush(fmass);
  fclose(fmass);
}

