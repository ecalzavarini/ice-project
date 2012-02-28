#include "common.h"


#ifdef FLUID_RHEOLOGY
void relaxation_time(){
  int x, y, pp, i, j;
  double u, v;
  double invtau, rho;
  double cu, u2;
  pop p_eq;
  double tau0;
  double S[2][2], gamma_dot,eps;
  double nu0 =  (tau1-0.5)/3.0;
  double nu00 = 0.5*(tau1-0.5)/3.0;  
  double lambda = 1.0;
  double nindex = 0.1;

  double B, tau_min;


  
  for (y=1; y<NY+1; y++){
    for (x=1; x<NX+1; x++){

      rho = m(p[IDX(y,x)]);
      u = vx(p[IDX(y,x)])/rho;
      v = vy(p[IDX(y,x)])/rho;
      u2 = u*u +  v*v;

      tau0 = tau[y][x];

      /* equilibrium distribution */
      for (pp=0; pp<9; pp++){
	cu = (cx[pp]*u + cy[pp]*v);
	p_eq.p[pp] = rho * wgt[pp] * (1.0 + 3.0*cu  + 4.5*cu*cu - 1.5*u2 );
      }

      S[0][0] = S[0][1] = S[1][0] = S[1][1] = 0.0;
      for (pp=0; pp<9; pp++){
	S[0][0] += cx[pp]*cx[pp]*(p[IDX(y,x)].p[pp] - p_eq.p[pp]);
	S[0][1] += cx[pp]*cy[pp]*(p[IDX(y,x)].p[pp] - p_eq.p[pp]);
	S[1][0] += cy[pp]*cx[pp]*(p[IDX(y,x)].p[pp] - p_eq.p[pp]);
	S[1][1] += cy[pp]*cy[pp]*(p[IDX(y,x)].p[pp] - p_eq.p[pp]);
      }

      //fprintf(stderr,"Sxx %g, Sxy %g, Syx %g, Syy %g\n", Sxx, Sxy, Syx, Syy);fflush(stderr);

      /* shear rate */
      eps = 0.0;
      for (i=0; i<3; i++) 
	for (j=0; j<3; j++){
	  eps += S[i][j]*S[i][j];
	}
      gamma_dot = (1.0/(2.0*cs2*rho*tau0))*sqrt(eps);

      /* fprintf(stderr,"gamma_dot %g, eps %e\n", gamma_dot, eps);fflush(stderr); */

      /* Carreau-Yasuda model , has parameters nu0, nu00 , lambda, nindex */
#ifdef FLUID_RHEOLOGY_CARREAU
      tau[IDX(y,x)]=((nu0 - nu00)/cs2) * pow((1.0 + pow(lambda*gamma_dot,2.0) ),(nindex-1.0)/2.0) + nu00/cs2 + 0.5;
#endif

      /* Power law model */
#ifdef FLUID_RHEOLOGY_POWER_LAW
      tau[IDX(y,x)]= (nu0/cs2) * pow(gamma_dot, nindex-1.0 ) + 0.5;
#endif

      /* if( tau[y][x] != tau0 ) fprintf(stderr,"tau %e tau0 %e\n", tau[y][x],tau0); */

    }
  }


#ifdef TEMP_REHOLOGY
  tau_min = 2./3.;
  B = deltaT;
  nu0 = (tau_min-0.5)*cs2;
  nu0 = nu0*exp(B/tt[IDX(y,x)]);

  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++){
      tau[IDX(y,x)] = nu0*exp(-B/tt[IDX(y,x)])/cs2 + 0.5;
    }
#endif

}
#endif
