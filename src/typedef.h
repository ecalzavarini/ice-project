#include "define.h"

//#define my_double long double 
#define my_double double

/* Useful structures */
typedef struct {
  my_double p[9];
} pop;

typedef struct {
  my_double val;
  my_double name;
} param;

typedef struct {
  int NX , NY;
#ifdef FLUID
  my_double tau;
  my_double nu;
#endif
#ifdef TEMPERATURE
  my_double tau_t;
  my_double kappa_t;
  my_double T_bot, T_top, T_ref;
  my_double deltaT;
#ifdef TEMPERATURE_BUOYANCY
  my_double beta_t;
  my_double beta2_t;
  my_double gravity_x;
  my_double gravity_y;
#endif
#ifdef TEMPERATURE_MELTING
  my_double T_solid;
  my_double specific_heat;
  my_double latent_heat;
  my_double liquidus_slope;  
#endif
#endif
#ifdef SALT
  my_double tau_s;
  my_double kappa_s;
  my_double S_bot, S_top, S_ref;
  my_double deltaS;
#ifdef SALT_BUOYANCY
  my_double beta_s;
  my_double beta2_s;
#ifndef TEMPERATURE_BUOYANCY
  my_double gravity_x;
  my_double gravity_y;
#endif
#endif
#endif
} prop;

typedef struct {
#ifdef TEMPERATURE
  my_double Prandtl;
#ifdef TEMPERATURE_BUOYANCY
  my_double Rayleigh_t;
#endif
#ifdef TEMPERATURE_MELTING
  my_double Stefan;
#endif
#endif
#ifdef SALT
  my_double Schmidt;
#endif
#ifdef SALT && TEMPERATURE
  my_double Lewis;
#endif
#ifdef SALT_BUOYANCY
  my_double Rayleigh_s;
#endif
} dimensionless_prop;

typedef struct {
  my_double vx;
  my_double vy;
} velocity;

typedef struct {
  my_double x;
  my_double y;
} vector;


typedef struct {
  my_double vx, vy, rho;
  my_double vx2, vy2;
#ifdef TEMPERATURE
  my_double nusselt, vyt;
  my_double dyt;
  my_double t , t2;
#endif
#ifdef SALT
  my_double s, s2 ;
#endif
#ifdef TEMPERATURE_MELTING
  my_double lf;
#endif
} diag;


/* WARNING vx means rho*vx */
#define vx(a) (a.p[1]+a.p[5]+a.p[8]-a.p[3]-a.p[6]-a.p[7])
#define vy(a) (a.p[2]+a.p[5]+a.p[6]-a.p[4]-a.p[7]-a.p[8])
#define  m(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8])
#define  t(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8])


#define IDX(j,i) ( (int)(j) * ( (NX) + 2 ) + (int)(i) )

#define IDXbuffer(plane,i) ( (int)(plane) * ( (NX) + 2 ) + (int)(i) )
