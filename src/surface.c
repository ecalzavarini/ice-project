#include "common_object.h"


#ifdef FLUID_SURFACE
void surface(int i){
  int x, y , pp;
  my_double  fl;
  my_double hs, hl , temp , Tl;
  my_double eps = 1.0;
  my_double Ts,Cp,Lf;

  /* store previous mass fraction */
  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++) {
      mmold[IDX(y,x)] = mm[IDX(y,x)]; 
    }

  /* compute mass fraction in each cell*/
  


  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++) {
      /* compute Entalphy */
      hh[IDX(y,x)] = Cp*tt[IDX(y,x)] + Lf*ll[IDX(y,x)];
      /* compute new fluid fraction */
      if(hh[IDX(y,x)] < hs) ll[IDX(y,x)]=0.0;
      else if(hh[IDX(y,x)] > hl) ll[IDX(y,x)]=1.0;
      else ll[IDX(y,x)]= (hh[IDX(y,x)] - hs)/(hl - hs);
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



}
#endif
