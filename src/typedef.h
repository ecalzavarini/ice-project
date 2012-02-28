#include "define.h"

/* Useful structures */
typedef struct {
  double p[9];
} pop;

typedef struct {
  double val;
  double name;
} param;

typedef struct {
  int NX , NY;
#ifdef FLUID
  double tau;
  double nu;
#endif
#ifdef TEMPERATURE
  double tau_t;
  double kappa_t;
  double T_bot, T_top, T_ref;
  double deltaT;
#ifdef TEMPERATURE_BUOYANCY
  double beta_t;
  double beta2_t;
  double gravity_x;
  double gravity_y;
#endif
#ifdef TEMPERATURE_MELTING
  double T_solid;
  double specific_heat;
  double latent_heat;
  double liquidus_slope;  
#endif
#endif
#ifdef SALT
  double tau_s;
  double kappa_s;
  double S_bot, S_top, S_ref;
  double deltaS;
#ifdef SALT_BUOYANCY
  double beta_s;
  double beta2_s;
#ifndef TEMPERATURE_BUOYANCY
  double gravity_x;
  double gravity_y;
#endif
#endif
#endif
} prop;

typedef struct {
#ifdef TEMPERATURE
  double Prandtl;
#ifdef TEMPERATURE_BUOYANCY
  double Rayleigh_t;
#ifdef TEMPERATURE_MELTING
  double Stefan;
#endif
#endif
#endif
#ifdef SALT
  double Schmidt;
#endif
#ifdef SALT && TEMPERATURE
  double Lewis;
#endif
#ifdef SALT_BUOYANCY
  double Rayleigh_s;
#endif
} dimensionless_prop;

typedef struct {
  double vx;
  double vy;
} velocity;

typedef struct {
  double x;
  double y;
} vector;


typedef struct {
  double vx, vy, rho;
  double vx2, vy2;
#ifdef TEMPERATURE
  double nusselt, vyt;
  double dyt;
  double t , t2;
#endif
#ifdef SALT
  double s, s2 ;
#endif
#ifdef TEMPERATURE_MELTING
  double lf;
#endif
} diag;


/* WARNING vx means rho*vx */
#define vx(a) (a.p[1]+a.p[5]+a.p[8]-a.p[3]-a.p[6]-a.p[7])
#define vy(a) (a.p[2]+a.p[5]+a.p[6]-a.p[4]-a.p[7]-a.p[8])
#define  m(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8])
#define  t(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8])


#define IDX(j,i) ( (int)(j) * ( (NX) + 2 ) + (int)(i) )

#define IDXbuffer(plane,i) ( (int)(plane) * ( (NX) + 2 ) + (int)(i) )
