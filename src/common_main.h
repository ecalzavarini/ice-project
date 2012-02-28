#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include "functions.h"

/* System size */
int NX , NY;

/* Array of properties*/
prop property;
dimensionless_prop dimensionless;
diag rbout;

#ifdef FLUID
/* Populations */
pop *p;
pop *buffer;
velocity *v,*vold;
double *dens;
vector *force;
#endif

/* LB speeds & weights */
double wgt[9] , cx[9] , cy[9]; 
double dirp[9] ,invp[9];
double cs , cs2 , cs4 , twocs2 , twocs4;
double invcs , invcs2 ,  invcs4 , invtwocs2 , invtwocs4;

/* time */
int itime;
int max_step;
int time_dump_field;
char OutDir[256];

#ifdef TEMPERATURE
/* Populations */
pop *g;
double *tt,*ttold;
#endif

#ifdef SALT
/* Populations */
pop *s;
double *ss,*ssold;
#endif

#ifdef TEMPERATURE_MELTING
double *ll,*llold, *hh;
#endif

#ifdef FLUID_RHEOLOGY
double *tau;
#endif














