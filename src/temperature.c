#include "common_object.h"

//extern pop equilibrium();

#ifdef TEMPERATURE
void bct() 
{
  int pp;
  int x, y;
  int xm, xp;
  my_double rho,rhot,ux,uy,cu,u2; 
  pop g_eq, g_eq_w;
  my_double effDT; /*effective DT */

  for (x=0; x<NX+2; x++) {
    for (pp=0; pp<9; pp++) {
      g[IDX(0,x)].p[pp] = 0.0;
      g[IDX(NY+1,x)].p[pp] = 0.0;
    }
  }

  for (y=0; y<NY+2; y++) {
    for (pp=0; pp<9; pp++) {
      g[IDX(y,0)].p[pp] = 0.0;
      g[IDX(y,NX+1)].p[pp] = 0.0;
    }
  }

#ifdef TEMPERATURE_BC_ISO_Y
  /* bounce back */
  for (x=1; x<NX+1; x++) {
    xm = x-1;
    xp = x+1;

    /* equilibrium distribution in y=1 */
    /*
#ifdef METHOD_FORCING_GUO
    rhot =  tt[IDX(1,x)];
    rho = dens[IDX(1,x)];
     ux = (v[IDX(1,x)].vx + 0.5*force[IDX(1,x)].x)/dens[IDX(1,x)];
     uy = (v[IDX(1,x)].vy + 0.5*force[IDX(1,x)].y)/dens[IDX(1,x)];
#else
    */
    //rho = m(p[IDX(1,x)]);
    //rhot =  tt[IDX(1,x)];
    rhot = t(g[IDX(1,x)]);
    ux = v[IDX(1,x)].vx;
    uy = v[IDX(1,x)].vy;
    /*
#endif
    */
    u2 = ux*ux+uy*uy;    
    for (pp=0; pp<9; pp++){
      cu = (cx[pp]*ux + cy[pp]*uy);
      g_eq.p[pp] = rhot * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*cu*cu - invtwocs2*u2 );
    } 
      

    /*     effDT = ( deltaT/2.0 - t(g_eq) )*2.0 +  t(g_eq); */
    effDT = ( (property.T_bot-property.T_ref) - t(g_eq) )*2.0 +  t(g_eq);

    for (pp=0; pp<9; pp++){
      cu = -(cx[pp]*ux + cy[pp]*uy);
      g_eq_w.p[pp] = effDT * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*cu*cu - invtwocs2*u2 );
    } 

    /* to be tested
    g_eq = equilibrium(g,1,x);      

    effDT = ( (property.T_bot-property.T_ref) - t(g_eq) )*2.0 +  t(g_eq);
    rhot = t(g[IDX(y,x)]);

    g_eq_w = equilibrium_opposite_velocity(g,1,x);

    for (pp=0; pp<9; pp++)  g_eq_w.p[pp] = (g_eq_w.p[pp]/rhot)*effDT; 
    */

    /* At bottom */
#ifdef TEMPERATURE_BC_ISO_Y_BBACK
#ifdef MIRROR
    /* mirror */
    /*
    g[0][x ].p[2] = g[1][x].p[4] - g_eq.p[4] + wgt[2]*deltaT/2.0; 
    g[0][xm].p[5] = g[1][x].p[8] - g_eq.p[8] + wgt[5]*deltaT/2.0; 
    g[0][xp].p[6] = g[1][x].p[7] - g_eq.p[7] + wgt[6]*deltaT/2.0;
    */
    /*
    g[0][x ].p[2] = g[1][x].p[4] - g_eq.p[4] + wgt[2]*effDT; 
    g[0][xm].p[5] = g[1][x].p[8] - g_eq.p[8] + wgt[5]*effDT; 
    g[0][xp].p[6] = g[1][x].p[7] - g_eq.p[7] + wgt[6]*effDT;
    */
    g[IDX(0,x )].p[2] = g[IDX(1,x)].p[4] - g_eq.p[4] + g_eq_w.p[2]; 
    g[IDX(0,xm)].p[5] = g[IDX(1,x)].p[8] - g_eq.p[8] + g_eq_w.p[5]; 
    g[IDX(0,xp)].p[6] = g[IDX(1,x)].p[7] - g_eq.p[7] + g_eq_w.p[6];
#else
    /* opposite */
    /*
    g[0][x ].p[2] = g[1][x].p[4] - g_eq.p[4] + wgt[2]*deltaT/2.0; 
    g[0][xm].p[5] = g[1][x].p[7] - g_eq.p[7] + wgt[5]*deltaT/2.0; 
    g[0][xp].p[6] = g[1][x].p[8] - g_eq.p[8] + wgt[6]*deltaT/2.0; 
    */
    /* opposite */
    /*
    g[0][x ].p[2] = g[1][x].p[4] - g_eq.p[4] + wgt[2]*effDT; 
    g[0][xm].p[5] = g[1][x].p[7] - g_eq.p[7] + wgt[5]*effDT; 
    g[0][xp].p[6] = g[1][x].p[8] - g_eq.p[8] + wgt[6]*effDT; 
    */
    g[IDX(0,x) ].p[2] = g[IDX(1,x)].p[4] - g_eq.p[4] + g_eq_w.p[2]; 
    g[IDX(0,xm)].p[5] = g[IDX(1,x)].p[7] - g_eq.p[7] + g_eq_w.p[5]; 
    g[IDX(0,xp)].p[6] = g[IDX(1,x)].p[8] - g_eq.p[8] + g_eq_w.p[6]; 
#endif

#else /* #ifdef TEMPERATURE_BC_ISO_Y_BBACK*/
     effDT =  ( (property.T_bot-property.T_ref) - rhot )*2.0 +  rhot;
for (pp=0; pp<9; pp++){
  // g[IDX(0,x) ].p[pp] = g_eq_w.p[pp];
     g[IDX(0,x)].p[pp] = effDT * wgt[pp];
 }
#endif

    /* equilibrium distribution in y=NY */
    /*
#ifdef METHOD_FORCING_GUO
     rho = dens[IDX(NY,x)];
     rhot = tt[IDX(NY,x)];
     ux = (v[IDX(NY,x)].vx + 0.5*force[IDX(NY,x)].x)/dens[IDX(NY,x)];
     uy = (v[IDX(NY,x)].vy + 0.5*force[IDX(NY,x)].y)/dens[IDX(NY,x)];
#else
    */
    //     rho = m(p[IDX(NY,x)]);
    //rhot = tt[IDX(NY,x)];
     rhot = t(g[IDX(NY,x)]);
      ux = v[IDX(NY,x)].vx;
      uy = v[IDX(NY,x)].vy;
      /*
#endif
      */
    u2 = ux*ux+uy*uy;    
    for (pp=0; pp<9; pp++){
      cu = (cx[pp]*ux + cy[pp]*uy);
      g_eq.p[pp] = rhot * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*cu*cu - invtwocs2*u2 );
    }


    /*   effDT =  (- deltaT/2.0 - t(g_eq) )*2.0 +  t(g_eq); */
    effDT =  ( (property.T_top-property.T_ref) - t(g_eq) )*2.0 +  t(g_eq);


    for (pp=0; pp<9; pp++){
      cu = -(cx[pp]*ux + cy[pp]*uy);
      g_eq_w.p[pp] = effDT * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*cu*cu - invtwocs2*u2 );
    } 
    /* Mistake
    g_eq = equilibrium(g,NY,x);   
    effDT =  ( (property.T_top-property.T_ref) - t(g_eq) )*2.0 +  t(g_eq);
    rhot = t(g[IDX(y,x)]);

    for (pp=0; pp<9; pp++)  g_eq_w.p[pp] = (g_eq.p[pp]/rhot)*effDT; 
    */

    /* at top */
#ifdef TEMPERATURE_BC_ISO_Y_BBACK
#ifdef MIRROR
    /* mirror */
    /*
    g[NY+1][x ].p[4] = g[NY][x].p[2] - g_eq.p[2] + wgt[4]*(-deltaT/2.0);  
    g[NY+1][xm].p[8] = g[NY][x].p[5] - g_eq.p[5] + wgt[8]*(-deltaT/2.0);  
    g[NY+1][xp].p[7] = g[NY][x].p[6] - g_eq.p[6] + wgt[7]*(-deltaT/2.0);

    g[NY+1][x ].p[4] = g[NY][x].p[2] - g_eq.p[2] + wgt[4]*effDT;  
    g[NY+1][xm].p[8] = g[NY][x].p[5] - g_eq.p[5] + wgt[8]*effDT;  
    g[NY+1][xp].p[7] = g[NY][x].p[6] - g_eq.p[6] + wgt[7]*effDT;
    */
    g[IDX(NY+1,x)].p[4] = g[IDX(NY,x)].p[2] - g_eq.p[2] + g_eq_w.p[4];  
    g[IDX(NY+1,xm)].p[8] = g[IDX(NY,x)].p[5] - g_eq.p[5] + g_eq_w.p[8];  
    g[IDX(NY+1,xp)].p[7] = g[IDX(NY,x)].p[6] - g_eq.p[6] + g_eq_w.p[7]; 
#else
    /* opposite */
    /*
    g[NY+1][x ].p[4] = g[NY][x].p[2] - g_eq.p[2] + wgt[4]*(-deltaT/2.0);  
    g[NY+1][xm].p[8] = g[NY][x].p[6] - g_eq.p[6] + wgt[8]*(-deltaT/2.0);  
    g[NY+1][xp].p[7] = g[NY][x].p[5] - g_eq.p[5] + wgt[7]*(-deltaT/2.0);  
    */
    /* opposite */
    /*
    g[NY+1][x ].p[4] = g[NY][x].p[2] - g_eq.p[2] + wgt[4]*effDT;  
    g[NY+1][xm].p[8] = g[NY][x].p[6] - g_eq.p[6] + wgt[8]*effDT;  
    g[NY+1][xp].p[7] = g[NY][x].p[5] - g_eq.p[5] + wgt[7]*effDT; 
    */
    g[IDX(NY+1,x) ].p[4] = g[IDX(NY,x)].p[2] - g_eq.p[2] + g_eq_w.p[4];  
    g[IDX(NY+1,xm)].p[8] = g[IDX(NY,x)].p[6] - g_eq.p[6] + g_eq_w.p[8];  
    g[IDX(NY+1,xp)].p[7] = g[IDX(NY,x)].p[5] - g_eq.p[5] + g_eq_w.p[7]; 
#endif    

#else /* #ifdef TEMPERATURE_BC_ISO_Y_BBACK */
 effDT =  ( (property.T_top-property.T_ref) - rhot )*2.0 +  rhot;
for (pp=0; pp<9; pp++){
  //g[IDX(NY+1,x)].p[pp] = g_eq_w.p[pp];
  g[IDX(NY+1,x)].p[pp] = effDT * wgt[pp];
 }
#endif

  }
#else  
  for (x=1; x<NX+1; x++) {
  /* periodic bc at top and bottom */
    g[IDX(0,x)] = g[IDX(1,x)];
    g[IDX(NY+1,x)] = g[IDX(NY,x)];
  }
#endif

  /* periodic bc at left and right */
  for (y=1; y<NY+1; y++) {
    g[IDX(y,   0)] = g[IDX(y,NX)];
    g[IDX(y,NX+1)] = g[IDX(y, 1)];
  }
}
#endif



#ifdef TEMPERATURE_BUOYANCY
#ifdef METHOD_FORCING_GUO
void buoyancy() 
{
  my_double ff1, ff2, vel, temp, coeff;
  int x, y;

  ff1 = property.beta_t*property.gravity_y;
  ff2 = property.beta2_t*property.gravity_y;

  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){

      temp = (tt[IDX(y,x)] - property.T_ref);
      force[IDX(y,x)].x += 0.0;
      force[IDX(y,x)].y += (ff1*temp + ff2*temp*temp);
      //if(itime==1)fprintf (stderr, "%g %g\n",  tt[IDX(y,x)], force[IDX(y,x)].y);
    }
}
#else
void buoyancy() 
{
  my_double ff1, ff2, vel, temp;
  int x, y, pp;

  ff1 = property.beta_t*property.gravity_y;
  ff2 = property.beta2_t*property.gravity_y;
  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){

      temp = 3.0* (tt[IDX(y,x)] - property.T_ref) * dens[IDX(y,x)];
      /*temp = 3.0*t(g[y][x]) * m(p[y][x]);*/ 
      for (pp=0; pp<9; pp++) p[IDX(y,x)].p[pp] += wgt[pp]*cy[pp]*(ff1*temp + ff2*temp*temp);
      /*
      p[y][x].p[6] += wgt[6]*ff_true*temp*cy[6];
      p[y][x].p[2] += wgt[2]*ff_true*temp*cy[2];
      p[y][x].p[5] += wgt[5]*ff_true*temp*cy[5];

      p[y][x].p[7] += wgt[7]*ff_true*temp*cy[7];
      p[y][x].p[4] += wgt[4]*ff_true*temp*cy[4];
      p[y][x].p[8] += wgt[8]*ff_true*temp*cy[8];
      */
    }
}
#endif
#endif


#ifdef TEMPERATURE_FORCING
void temperature_forcing() 
{
  int x , y, pp; 
  my_double internal_heat = 5.333333e-4;

  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){
      for (pp=0; pp<9; pp++) g[IDX(y,x)].p[pp] += wgt[pp]*internal_heat;
    }

}
#endif
