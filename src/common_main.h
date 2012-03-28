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
my_double *dens;
vector *force;
#endif

/* LB speeds & weights */
my_double wgt[9] , cx[9] , cy[9]; 
my_double dirp[9] ,invp[9];
my_double cs , cs2 , cs4 , twocs2 , twocs4;
my_double invcs , invcs2 ,  invcs4 , invtwocs2 , invtwocs4;

/* time */
int itime;
int max_step;
int time_dump_field;
char OutDir[256];

#ifdef TEMPERATURE
/* Populations */
pop *g;
my_double *tt,*ttold;
#endif

#ifdef SALT
/* Populations */
pop *s;
my_double *ss,*ssold;
#endif

#ifdef TEMPERATURE_MELTING
my_double *ll,*llold, *hh;
#endif

#ifdef FLUID_RHEOLOGY
my_double *tau;
#endif














