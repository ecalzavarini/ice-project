#include "common.h"


#ifdef FLUID_RHEOLOGY_POLYMER
void compute_fluid_gradient_tensor(){
  int x, y, pp, i, j;
  pop p_eq;
  tensor S;

  for (y=1; y<NY+1; y++){
    for (x=1; x<NX+1; x++){

      /* equilibrium distribution */
      p_eq=equilibrium(p,y,x);

      S.xx = S.xy = Sy.x = S.yy = 0.0;
      for (pp=0; pp<9; pp++){
	S.xx += cx[pp]*cx[pp]*(p[IDX(y,x)].p[pp] - p_eq.p[pp]);
	S.xy += cx[pp]*cy[pp]*(p[IDX(y,x)].p[pp] - p_eq.p[pp]);
	S.yx += cy[pp]*cx[pp]*(p[IDX(y,x)].p[pp] - p_eq.p[pp]);
	S.yy += cy[pp]*cy[pp]*(p[IDX(y,x)].p[pp] - p_eq.p[pp]);
      }
      gradv[IDX(y,x)]=S;

      //fprintf(stderr,"S.xx %g, S.xy %g, S.yx %g, S.yy %g\n", S.xx, S.xy, S.yx, S.yy);fflush(stderr);
    }
  }
}
#endif


tensor mat_multiply(tensor A, tensor B){
  tensor C;
  C.xx = A.xx*B.xx+ A.xy*B.yx;
  C.xy = A.xx*B.yx+ A.xy*B.yy;
  C.yx = A.yx*B.xx+ A.yy*B.yx;
  C.yy = A.yx*B.yx+ A.yy*B.yy;
  return C;
}

tensor mat_transpose(tensor A){
  tensor B;
  C.xx = B.xx;
  C.xy = B.xy;
  C.yx = B.yx;
  C.yy = B.yy;
  return B;
}

 

void compute_polymer_gradient_tensor(){
/* Compute first derivatives of the conformation matrix - finite difference*/
  int xp,yp,xm,ym;
  velocity vel;
  tensor dx_conf,dy_conf;
  for (y=1; y<NY+1; y++){
    for (x=1; x<NX+1; x++){      
      xp=x+1; xm=x-1;
      yp=y+1; ym=y-1;

      vel=v[IDX(y,x)];

      dx_conf.xx = 0.5*(conf[IDX(y,xp)].xx-conf[IDX(y,xm)].xx);
      dx_conf.xy = 0.5*(conf[IDX(y,xp)].xy-conf[IDX(y,xm)].xy);
      dx_conf.yy = 0.5*(conf[IDX(y,xp)].yy-conf[IDX(y,xm)].yy);
      dx_conf.yx = dx_conf.xy;

      dy_conf.xx = 0.5*(conf[IDX(yp,x)].xx-conf[IDX(ym,x)].xx);
      dy_conf.xy = 0.5*(conf[IDX(yp,x)].xy-conf[IDX(ym,x)].xy);
      dy_conf.yy = 0.5*(conf[IDX(yp,x)].yy-conf[IDX(ym,x)].yy);
      dy_conf.yx = dx_conf.xy;

      rhs_conf[IDX(y,x)].xx -= (vel.vx*dx_conf.xx + vel.vy*dx_conf.xx);
      rhs_conf[IDX(y,x)].xy -= (vel.vx*dx_conf.xy + vel.vy*dx_conf.xy);
      rhs_conf[IDX(y,x)].yx -= (vel.vx*dx_conf.yx + vel.vy*dx_conf.yx);
      rhs_conf[IDX(y,x)].yy -= (vel.vx*dx_conf.yy + vel.vy*dx_conf.yy);

    }
  }
}


void compute_polymer_laplacian(){
/* Compute laplacian of conformation matrix - finite difference*/
int xp,yp,xm,ym;
 tensor dxx_conf, dyy_conf;
  for (y=1; y<NY+1; y++){
    for (x=1; x<NX+1; x++){      
      xp=x+1; xm=x-1;
      yp=y+1; ym=y-1;

      dxx_conf.xx = (conf[IDX(y,xp)].xx-2.*conf[IDX(y,x)].xx+conf[IDX(y,xm)].xx);
      dxx_conf.xy = (conf[IDX(y,xp)].xy-2.*conf[IDX(y,x)].xy+conf[IDX(y,xm)].xy);
      dxx_conf.yy = (conf[IDX(y,xp)].yy-2.*conf[IDX(y,x)].yy+conf[IDX(y,xm)].yy);
      dxx_conf.yx = dxx_conf.xy;

      dyy_conf.xx = (conf[IDX(yp,x)].xx-2.*conf[IDX(y,x)].xx+conf[IDX(ym,x)].xx);
      dyy_conf.xy = (conf[IDX(yp,x)].xy-2.*conf[IDX(y,x)].xy+conf[IDX(ym,x)].xy);
      dyy_conf.yy = (conf[IDX(yp,x)].yy-2.*conf[IDX(y,x)].yy+conf[IDX(ym,x)].yy);
      dyy_conf.yx = dyy_conf.xy;

      rhs_conf[IDX(y,x)].xx += dxx_conf.xx + dyy_conf.xx; 
      rhs_conf[IDX(y,x)].xy += dxx_conf.xy + dyy_conf.xy;
      rhs_conf[IDX(y,x)].yy += dxx_conf.yy + dyy_conf.yy;
      rhs_conf[IDX(y,x)].yx += dxx_conf.yx + dyy_conf.yx;
    }
  }
}


void compute_elastic_term(){
  double two_invtau_polymer = 2.0/tau_polymer;

  for (y=1; y<NY+1; y++){
    for (x=1; x<NX+1; x++){            
   /* elastic term */
   rhs_conf[IDX(y,x)].xx += two_invtau_polymer*(conf[IDX(y,x)].xx - 1.0);
   rhs_conf[IDX(y,x)].yy += two_invtau_polymer*(conf[IDX(y,x)].yy - 1.0);
   rhs_conf[IDX(y,x)].xy += two_invtau_polymer*conf[IDX(y,x)].xy; 
   rhs_conf[IDX(y,x)].yx += two_invtau_polymer*conf[IDX(y,x)].yx;
    }
  }
}


/* Assemble right-hand-side */

/* compute: \nabla u */
   compute_fluid_gradient_tensor();

/* compute: - (u \cdot \nabla)\sigma */
   compute_polymer_gradient_tensor();

/* compute:  \sigma* \nabla u  */ 
   rhs_conf = mat_multiply(conf,gradv);

/* compute:  (\nabla u)^T * \sigma  */
     gradvT = mat_transpose(gradv);
   rhs_conf = mat_multiply(gradvT,conf);

/* compute:  (2 / \tau_polymer)*(\sigma - Identity)  */
   compute_elastic_term();

/* compute:  \kappa_polymer * \Delta \sigma  */
   compute_polymer_laplacian();
