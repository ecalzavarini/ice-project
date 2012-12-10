#include "common.h"




#ifdef FLUID_RHEOLOGY_POLYMER
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
#endif
 

#ifdef FLUID_RHEOLOGY_POLYMER
void compute_polymer_gradient_tensor(){
/* Compute first derivatives of the conformation matrix - finite difference*/
  int xp,yp,xm,ym;
  velocity vel;
  tensor dx_conf,dy_conf;
  double fac;
  for (y=1; y<NY+1; y++){
    for (x=1; x<NX+1; x++){      

      vel=v[IDX(y,x)];

      if(x>1 && x<NX+1){
	xp=x+1; xm=x-1;
	fac=0.5;
      }else{
      if(x==1){
	xp=x+1; xm=x;
	fac=1.0;
      }
      if(x==NX+1){
	xp=x; xm=x-1;
	fac=1.0;
      }      
      }

      dx_conf.xx = fac*(conf[IDX(y,xp)].xx-conf[IDX(y,xm)].xx);
      dx_conf.xy = fac*(conf[IDX(y,xp)].xy-conf[IDX(y,xm)].xy);
      dx_conf.yy = fac*(conf[IDX(y,xp)].yy-conf[IDX(y,xm)].yy);
      dx_conf.yx = dx_conf.xy;

      if(y>1 && y<NY+1){
	yp=y+1; ym=y-1;
	fac=0.5;
      }else{
      if(y==1){
	yp=y+1; ym=y;
	fac=1.0;
      }
      if(y==NY+1){
	yp=y; ym=y-1;
	fac=1.0;
      }        
      }

      dy_conf.xx = fac*(conf[IDX(yp,x)].xx-conf[IDX(ym,x)].xx);
      dy_conf.xy = fac*(conf[IDX(yp,x)].xy-conf[IDX(ym,x)].xy);
      dy_conf.yy = fac*(conf[IDX(yp,x)].yy-conf[IDX(ym,x)].yy);
      dy_conf.yx = dy_conf.xy;

      rhs_conf[IDX(y,x)].xx -= (vel.vx*dx_conf.xx + vel.vy*dx_conf.xx);
      rhs_conf[IDX(y,x)].xy -= (vel.vx*dx_conf.xy + vel.vy*dx_conf.xy);
      rhs_conf[IDX(y,x)].yx -= (vel.vx*dx_conf.yx + vel.vy*dx_conf.yx);
      rhs_conf[IDX(y,x)].yy -= (vel.vx*dx_conf.yy + vel.vy*dx_conf.yy);

    }
  }
}
#endif


#ifdef FLUID_RHEOLOGY_POLYMER
void compute_polymer_laplacian(){
/* Compute laplacian of conformation matrix - finite difference*/
  int xp,yp,xm,ym,x0,y0;
 tensor dxx_conf, dyy_conf;
  for (y=1; y<NY+1; y++){
    for (x=1; x<NX+1; x++){      

      /* x stencil */
      if(x>1 && x<NX+1){
	xp=x+1; x0=x; xm=x-1;
      }else{
      if(x==1){
	xp=x+2; x0=x+1; xm=x;
      }
      if(x==NX+1){
	xp=x; x0=x-1; xm=x-2;
      }      
      }
 
      dxx_conf.xx = (conf[IDX(y,xp)].xx-2.*conf[IDX(y,x0)].xx+conf[IDX(y,xm)].xx);
      dxx_conf.xy = (conf[IDX(y,xp)].xy-2.*conf[IDX(y,x0)].xy+conf[IDX(y,xm)].xy);
      dxx_conf.yy = (conf[IDX(y,xp)].yy-2.*conf[IDX(y,x0)].yy+conf[IDX(y,xm)].yy);
      dxx_conf.yx = dxx_conf.xy;

      /* y stencil */
      if(y>1 && y<NY+1){
	yp=y+1; y0=y; ym=y-1;
      }else{
      if(y==1){
	yp=y+2; y0=y+1; ym=y;
      }
      if(y==NY+1){
	yp=y; y0=y-1; ym=y-2;
      }      
      }

      dyy_conf.xx = (conf[IDX(yp,x)].xx-2.*conf[IDX(y0,x)].xx+conf[IDX(ym,x)].xx);
      dyy_conf.xy = (conf[IDX(yp,x)].xy-2.*conf[IDX(y0,x)].xy+conf[IDX(ym,x)].xy);
      dyy_conf.yy = (conf[IDX(yp,x)].yy-2.*conf[IDX(y0,x)].yy+conf[IDX(ym,x)].yy);
      dyy_conf.yx = dyy_conf.xy;

      rhs_conf[IDX(y,x)].xx += dxx_conf.xx + dyy_conf.xx; 
      rhs_conf[IDX(y,x)].xy += dxx_conf.xy + dyy_conf.xy;
      rhs_conf[IDX(y,x)].yy += dxx_conf.yy + dyy_conf.yy;
      rhs_conf[IDX(y,x)].yx += dxx_conf.yx + dyy_conf.yx;
    }
  }
}
#endif

#ifdef FLUID_RHEOLOGY_POLYMER
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
#endif

#ifdef FLUID_RHEOLOGY_POLYMER
void compute_time_step_euler(){
  for (y=1; y<NY+1; y++){
    for (x=1; x<NX+1; x++){            
   /*euler time step*/
     conf[IDX(y,x)].xx += rhs_conf[IDX(y,x)].xx; 
     conf[IDX(y,x)].yy += rhs_conf[IDX(y,x)].yy; 
     conf[IDX(y,x)].xy += rhs_conf[IDX(y,x)].xy; 
     conf[IDX(y,x)].yx += rhs_conf[IDX(y,x)].yx; 
    }
  }
}
#endif


#ifdef FLUID_RHEOLOGY_POLYMER
#ifdef FLUID_RHEOLOGY_POLYMER_FEEDBACK
void compute_polymer_extra_stress(){
/* Compute first derivatives of the conformation matrix - finite difference*/
  int xp,yp,xm,ym;
  velocity vel;
  tensor dx_conf,dy_conf;
  double fac;
  for (y=1; y<NY+1; y++){
    for (x=1; x<NX+1; x++){      

      vel=v[IDX(y,x)];

      /* x stencil */
      if(x>1 && x<NX+1){
	xp=x+1; xm=x-1;
	fac=0.5;
      }else{
      if(x==1){
	xp=x+1; xm=x;
	fac=1.0;
      }
      if(x==NX+1){
	xp=x; xm=x-1;
	fac=1.0;
      }      
      }

     /* y stencil */
      if(y>1 && y<NY+1){
	yp=y+1; ym=y-1;
	fac=0.5;
      }else{
      if(y==1){
	yp=y+1; ym=y;
	fac=1.0;
      }
      if(y==NY+1){
	yp=y; ym=y-1;
	fac=1.0;
      }        
      }

      dx_conf.xx = fac*(conf[IDX(y,xp)].xx-conf[IDX(y,xm)].xx);
      dy_conf.yx = fac*(conf[IDX(yp,x)].yx-conf[IDX(ym,x)].yx);

      dx_conf.xy = fac*(conf[IDX(y,xp)].xy-conf[IDX(y,xm)].xy);
      dy_conf.yy = fac*(conf[IDX(yp,x)].yy-conf[IDX(ym,x)].yy);

      force[IDX(y,x)].x += dx_conf.xx + dy_conf.yx;
      force[IDX(y,x)].y += dx_conf.xy + dy_conf.yy;

    }
  }
}
#endif
#endif

#ifdef FLUID_RHEOLOGY_POLYMER
void  polymer_conformation_equation(){
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

/* compute:  euler first order */
   compute_time_step_euler();

#ifdef FLUID_RHEOLOGY_POLYMER_FEEDBACK
/* compute the divergence  (\nabla * \sigma)  */
   compute_polymer_extra_stress();
#endif
}
#endif
