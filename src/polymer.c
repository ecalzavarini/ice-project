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


tensor mat_mul(tensor A, tensor B){
  C.xx = A.xx*B.xx+ A.xy*B.yx;
  C.xy = A.xx*B.yx+ A.xy*B.yy;
  C.yx = A.yx*B.xx+ A.yy*B.yx;
  C.yy = A.yx*B.yx+ A.yy*B.yy;
  return C;
}


/* compute left hand side */


mat_mul(gradv,conf)
mat_mul(conf,gradv)

