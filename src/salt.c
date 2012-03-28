#include "common_object.h"

#ifdef SALT
void bcs() 
{
  int pp;
  int x, y;
  int xm, xp;
  my_double rho,rhos,ux,uy,cu,u2; 
  pop s_eq, s_eq_w;
  my_double effDS; /*effective DS */

  for (x=0; x<NX+2; x++) {
    for (pp=0; pp<9; pp++) {
      s[IDX(0,x)].p[pp] = 0.0;
      s[IDX(NY+1,x)].p[pp] = 0.0;
    }
  }

  for (y=0; y<NY+2; y++) {
    for (pp=0; pp<9; pp++) {
      s[IDX(y,0)].p[pp] = 0.0;
      s[IDX(y,NX+1)].p[pp] = 0.0;
    }
  }


#ifdef SALT_BC_ISO_Y
  /* bounce back */
  for (x=1; x<NX+1; x++) {
    xm = x-1;
    xp = x+1;

    /* equilibrium distribution in y=1 */
    //rho = m(p[IDX(1,x)]);
    rhos = t(s[IDX(1,x)]);
    ux = v[IDX(1,x)].vx;
    uy = v[IDX(1,x)].vy;
    u2 = ux*ux+uy*uy;    
    for (pp=0; pp<9; pp++){
      cu = (cx[pp]*ux + cy[pp]*uy);
      s_eq.p[pp] = rhos * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*cu*cu - invtwocs2*u2 );
    } 
      

    /*     effDS = ( deltaS/2.0 - t(s_eq) )*2.0 +  t(s_eq);*/
    effDS = ( (property.S_bot-property.S_ref) - t(s_eq) )*2.0 +  t(s_eq);

    for (pp=0; pp<9; pp++){
      cu = -(cx[pp]*ux + cy[pp]*uy);
      s_eq_w.p[pp] = effDS * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*cu*cu - invtwocs2*u2 );
    } 


    /* At bottom */
#ifdef MIRROR
    /* mirror */
    s[IDX(0,x )].p[2] = s[IDX(1,x)].p[4] - s_eq.p[4] + s_eq_w.p[2]; 
    s[IDX(0,xm)].p[5] = s[IDX(1,x)].p[7] - s_eq.p[7] + s_eq_w.p[5]; 
    s[IDX(0,xp)].p[6] = s[IDX(1,x)].p[8] - s_eq.p[8] + s_eq_w.p[6];
#else
    /* opposite */
    s[IDX(0,x) ].p[2] = s[IDX(1,x)].p[4] - s_eq.p[4] + s_eq_w.p[2]; 
    s[IDX(0,xm)].p[5] = s[IDX(1,x)].p[7] - s_eq.p[7] + s_eq_w.p[5]; 
    s[IDX(0,xp)].p[6] = s[IDX(1,x)].p[8] - s_eq.p[8] + s_eq_w.p[6]; 
#endif

    /* equilibrium distribution in y=NY */
    //rho = m(p[IDX(NY,x)]);
    rhos = t(s[IDX(NY,x)]);
    ux = v[IDX(NY,x)].vx;
    uy = v[IDX(NY,x)].vy;
    u2 = ux*ux+uy*uy;    
    for (pp=0; pp<9; pp++){
      cu = (cx[pp]*ux + cy[pp]*uy);
      s_eq.p[pp] = rhos * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*cu*cu - invtwocs2*u2 );
    }


    /*     effDS =  (- deltaS/2.0 - t(s_eq) )*2.0 +  t(s_eq); */
    effDS =  ( (property.S_top-property.S_ref) - t(s_eq) )*2.0 +  t(s_eq);

    for (pp=0; pp<9; pp++){
      cu = -(cx[pp]*ux + cy[pp]*uy);
      s_eq_w.p[pp] = effDS * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*cu*cu - invtwocs2*u2 );
    } 

    /* at top */
#ifdef MIRROR
    /* mirror */
    g[IDX(NY+1,x) ].p[4] = g[IDX(NY,x)].p[2] - g_eq.p[2] + g_eq_w.p[4];  
    g[IDX(NY+1,xm)].p[8] = g[IDX(NY,x)].p[6] - g_eq.p[6] + g_eq_w.p[8];  
    g[IDX(NY+1,xp)].p[7] = g[IDX(NY,x)].p[5] - g_eq.p[5] + g_eq_w.p[7]; 
#else
    /* opposite */
    s[IDX(NY+1,x) ].p[4] = s[IDX(NY,x)].p[2] - s_eq.p[2] + s_eq_w.p[4];  
    s[IDX(NY+1,xm)].p[8] = s[IDX(NY,x)].p[6] - s_eq.p[6] + s_eq_w.p[8];  
    s[IDX(NY+1,xp)].p[7] = s[IDX(NY,x)].p[5] - s_eq.p[5] + s_eq_w.p[7]; 
#endif    
  }
#else
  for (x=1; x<NX+1; x++) {
  /* periodic bc at top and bottom */
    s[IDX(0,x)] = s[IDX(1,x)];
    s[IDX(NY+1,x)] = s[IDX(NY,x)];
  }
#endif

  /* periodic bc at left and right */
  for (y=1; y<NY+1; y++) {
    s[IDX(y,0)] = s[IDX(y,NX)];
    s[IDX(y,NX+1)] = s[IDX(y,1)];
  }
}
#endif




#ifdef SALT_BUOYANCY
#ifdef METHOD_FORCING_GUO
void buoyancys() 
{
  my_double ff_true, temp;
  int x, y;

  ff_true = property.beta_s*property.gravity_y;

  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){

      temp = (ss[IDX(y,x)] - property.S_ref);
      force[IDX(y,x)].x += 0.0;
      force[IDX(y,x)].y += ff_true*temp;
    }
}
#else
void buoyancys() 
{
  my_double ff_true, temp;
  int x, y, pp;

  ff_true = property.beta_s*property.gravity_y;
  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){
      temp = 3.0* ss[IDX(y,x)] * dens[IDX(y,x)];
      for (pp=0; pp<9; pp++) p[IDX(y,x)].p[pp] += wgt[pp]*ff_true*temp*cy[pp];
    }
}
#endif
#endif
