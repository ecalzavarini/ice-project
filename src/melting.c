#include "common_object.h"


#ifdef TEMPERATURE_MELTING
void melting(int i){
  int x, y , pp;
  double  fl;
  double hs, hl , temp , Tl;
  double eps = 1.0;
  double Ts,Cp,Lf;

  Ts = (double) property.T_solid;
  Cp = (double) property.specific_heat;
  Lf = (double) property.latent_heat;
  /* solidification temperature */
  hs = Cp*Ts;
  hl = hs+Lf;
  Tl = hl/Cp;

  if(i==0) fprintf(stderr,"Ts = %g \n Tl = %g\n",Ts, Tl); 

  /* store previous fluid fraction */
  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++) {
      llold[IDX(y,x)] = ll[IDX(y,x)]; 
    }

  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++) {
#ifdef SALT
      /* subtract the freezing point depression */
      Ts = property.liquidus_slope * ss[IDX(y,x)];
      hs = Cp*Ts;
      hl = hs+Lf;
#endif
      /* compute Entalphy */
      hh[IDX(y,x)] = Cp*tt[IDX(y,x)] + Lf*ll[IDX(y,x)];
      /* compute new fluid fraction */
      if(hh[IDX(y,x)] < hs) ll[IDX(y,x)]=0.0;
      else if(hh[IDX(y,x)] > hl) ll[IDX(y,x)]=1.0;
      else ll[IDX(y,x)]= (hh[IDX(y,x)] - hs)/(hl - hs);
    }

  /* add melting term to the temperature field */
  temp = (Lf/Cp);
  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){
      for (pp=0; pp<9; pp++) g[IDX(y,x)].p[pp] -= wgt[pp]*temp*(ll[IDX(y,x)]-llold[IDX(y,x)]);
    }

#ifdef METHOD_FORCING_GUO
  /* add drag term to the force structure */
  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){
      fl = ll[IDX(y,x)]; 
      temp = (1.0 - fl*fl)/(eps + fl*fl*fl);
      force[IDX(y,x)].x -= temp*v[IDX(y,x)].vx;
      force[IDX(y,x)].y -= temp*v[IDX(y,x)].vy;  
    }
#else
  /* add drag term to the velocity field */
  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){
      fl = ll[IDX(y,x)]; 
      temp = (1.0 - fl*fl)/(eps + fl*fl*fl);
      for (pp=0; pp<9; pp++) p[IDX(y,x)].p[pp] += -wgt[pp]*temp*(cx[pp]*v[IDX(y,x)].vx + cy[pp]*v[IDX(y,x)].vy);
    }
#endif




#ifdef THERMIC_SOLID
  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){
      fl = ll[IDX(y,x)];
      tauT[IDX(y,x)] = fl*tau2 + (1.0-fl)*tau4; 
    }
#endif

}
#endif
