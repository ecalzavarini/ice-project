#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include "define.h"
#include "typedef.h"

/* System size */
extern int NX , NY;

/* Array of properties*/
extern prop property;
extern dimensionless_prop dimensionless;
extern diag rbout;

#ifdef FLUID
/* Populations */
extern pop *p;
extern pop *buffer;
extern velocity *v,*vold;
extern double *dens;
extern vector *force;
#endif

/* LB speeds & weights */
extern double wgt[9] , cx[9] , cy[9]; 
extern double dirp[9] ,invp[9];
extern double cs , cs2 , cs4 , twocs2 , twocs4;
extern double invcs , invcs2 ,  invcs4 , invtwocs2 , invtwocs4;


/* time */
extern int itime;
extern int max_step;
extern int time_dump_field;
extern char OutDir[256];

#ifdef TEMPERATURE
/* Populations */
extern pop *g;
extern double *tt,*ttold;
#endif

#ifdef SALT
/* Populations */
extern pop *s;
extern double *ss,*ssold;
#endif

#ifdef TEMPERATURE_MELTING
extern double *ll,*llold, *hh;
#endif

#ifdef FLUID_RHEOLOGY
extern double *tau;
#endif














