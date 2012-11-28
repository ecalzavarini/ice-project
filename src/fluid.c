#include "common_object.h"

void bc() 
{
  int pp;
  int x, y;
  int xm, xp;
  int ym, yp;

  /* cornice superiore e inferiore -> set to zero */
  for (x=0; x<NX+2; x++) {
    for (pp=0; pp<9; pp++) {
      p[IDX(0,x)].p[pp] = 0.0;
      p[IDX(NY+1,x)].p[pp] = 0.0;
    }
  }

  /* cornice destra e sinistra -> set to zero */
  for (y=0; y<NY+2; y++) {
    for (pp=0; pp<9; pp++) {
      p[IDX(y,0)].p[pp] = 0.0;
      p[IDX(y,NX+1)].p[pp] = 0.0;
    }
  }


 /*********************************************                                                               
 *    y                                       *
 *    ^                                       *
 *    |                                       *
 *    |                                       *
 *    |                                       *
 *    +----------> x                          *
 *                                            *                                           
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
#ifdef FLUID_BC_WALL_Y
  for (x=1; x<NX+1; x++) {
    xm = x-1;
    xp = x+1;

#ifdef FLUID_BC_WALL_YM_SLIP
    /* free-slip bottom */
    p[IDX(0,x)].p[2] =  p[IDX(1,x)].p[4]; 
    p[IDX(0,xm)].p[5] = p[IDX(1,x)].p[6]; 
    p[IDX(0,xp)].p[6] = p[IDX(1,x)].p[5];
#else
    /* (bounce back) no-slip at bottom */ 
    p[IDX(0,x)].p[2] = p[IDX(1,x)].p[4]; 
    p[IDX(0,xm)].p[5] = p[IDX(1,x)].p[7]; 
    p[IDX(0,xp)].p[6] = p[IDX(1,x)].p[8]; 
#endif

#ifdef FLUID_BC_WALL_YP_SLIP
    /* free-slip at at top */
    p[IDX(NY+1,x) ].p[4] = p[IDX(NY,x)].p[2]; 
    p[IDX(NY+1,xm)].p[8] = p[IDX(NY,x)].p[7]; 
    p[IDX(NY+1,xp)].p[7] = p[IDX(NY,x)].p[8]; 
#else
    /* no-slip at at top */    
    p[IDX(NY+1,x) ].p[4] = p[IDX(NY,x)].p[2]; 
    p[IDX(NY+1,xp)].p[7] = p[IDX(NY,x)].p[5];  
    p[IDX(NY+1,xm)].p[8] = p[IDX(NY,x)].p[6];   
#endif
  }
#else
  for (x=1; x<NX+1; x++) {
/* periodic bc at top and bottom */
    p[IDX(0,x)] = p[IDX(NY,x)];
    p[IDX(NY+1,x)] = p[IDX(1,x)];
  }
#endif

 /*********************************************                                                               
 *    y                                       *
 *    ^                                       *
 *    |                                       *
 *    |                                       *
 *    |                                       *
 *    +----------> x                          *
 *                                            *                                           
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
#ifdef FLUID_BC_WALL_X
  for (y=1; y<NY+1; y++) {
    ym = y-1;
    yp = y+1;

#ifdef FLUID_BC_WALL_XM_SLIP
    /* free-slip left */
    p[IDX(y,0)].p[1] =  p[IDX(y,1)].p[3]; 
    p[IDX(ym,0)].p[8] = p[IDX(y,1)].p[7]; 
    p[IDX(yp,0)].p[5] = p[IDX(y,1)].p[6];
#else
    /* (bounce back) no-slip at left */
    p[IDX(y,0)].p[1] = p[IDX(y,1)].p[3]; 
    p[IDX(ym,0)].p[5] = p[IDX(y,1)].p[7]; 
    p[IDX(yp,0)].p[8] = p[IDX(y,1)].p[6]; 
#endif

#ifdef FLUID_BC_WALL_XP_SLIP
    /* free-slip at right */
    p[IDX(y,NX+1) ].p[3] = p[IDX(y,NX)].p[1]; 
    p[IDX(ym,NX+1)].p[7] = p[IDX(y,NX)].p[8]; 
    p[IDX(yp,NX+1)].p[6] = p[IDX(y,NX)].p[5]; 
#else
    /* no-slip at right */
    p[IDX(y,NX+1) ].p[3] = p[IDX(y,NX)].p[1]; 
    p[IDX(ym,NX+1)].p[6] = p[IDX(y,NX)].p[8]; 
    p[IDX(yp,NX+1)].p[7] = p[IDX(y,NX)].p[5]; 
#endif
  }
#else  
  for (y=1; y<NY+1; y++) {
  /* periodic bc at left and right */
    p[IDX(y,0)] = p[IDX(y,NX)];
    p[IDX(y,NX+1)] = p[IDX(y,1)];
  }
#endif

#ifdef FLUID_BC_FLOW_Y
  for (x=1; x<NX+1; x++) {
    xm = x-1;
    xp = x+1;

#ifdef FLUID_BC_FLOW_YM
    /* inflow at bottom */ 
    p[IDX(0,x)].p[2]  += ff_true;  p[IDX(1,x)].p[4] -= ff_true;
    p[IDX(0,xm)].p[5] += ff_true;  p[IDX(1,x)].p[7] -= ff_true; 
    p[IDX(0,xp)].p[6] += ff_true;  p[IDX(1,x)].p[8] -= ff_true; 
#endif

#ifdef FLUID_BC_FLOW_YP
    /* no-slip at at top */    
    p[IDX(NY+1,x) ].p[4] -= ff_true; p[IDX(NY,x)].p[2] += ff_true;; 
    p[IDX(NY+1,xp)].p[7] -= ff_true; p[IDX(NY,x)].p[5] += ff_true;;  
    p[IDX(NY+1,xm)].p[8] -= ff_true; p[IDX(NY,x)].p[6] += ff_true;;   
#endif
  }
#endif

}



#ifdef FLUID_FORCING_POISEUILLE
void poiseuille_forc(){
  double ff_true = property.gradP/6.0; //= 1.e-6;
  int x, y;

#ifdef METHOD_FORCING_GUO
  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++) {
      force[IDX(y,x)].x += ff_true;
    }
#else
  ff_true = (1./3.)*property.gradP/6.0;
  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++) {
      p[IDX(y,x)].p[1] += ff_true;
      p[IDX(y,x)].p[5] += ff_true;
      p[IDX(y,x)].p[8] += ff_true;

      p[IDX(y,x)].p[3] -= ff_true;
      p[IDX(y,x)].p[6] -= ff_true;
      p[IDX(y,x)].p[7] -= ff_true;
    }
#endif /* not method of forging guo */
}
#endif


#ifdef FLUID_FORCING_SPONGE
void sponge_forc(){
  my_double temp;
  int x, y;
  my_double vx0 = 0.0;
  my_double vy0 = 0.1;

  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){
#ifdef METHOD_FORCING_GUO
  /* add drag term to the force structure */
      temp = 
	force[IDX(y,x)].x -= temp*(v[IDX(y,x)].vx - vx0);
      force[IDX(y,x)].y -= temp*(v[IDX(y,x)].vy - vy0); 
#else
  /* add drag term to the velocity field */
      for (pp=0; pp<9; pp++) p[IDX(y,x)].p[pp] += -wgt[pp]*temp*(cx[pp]*(v[IDX(y,x)].vx - vx0) + cy[pp]*(v[IDX(y,x)].vy -vy0));
#endif
    }

}
#endif
