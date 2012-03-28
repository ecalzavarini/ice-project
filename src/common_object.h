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
extern my_double *dens;
extern vector *force;
#endif

/* LB speeds & weights */
extern my_double wgt[9] , cx[9] , cy[9]; 
extern my_double dirp[9] ,invp[9];
extern my_double cs , cs2 , cs4 , twocs2 , twocs4;
extern my_double invcs , invcs2 ,  invcs4 , invtwocs2 , invtwocs4;


/* time */
extern int itime;
extern int max_step;
extern int time_dump_field;
extern char OutDir[256];

#ifdef TEMPERATURE
/* Populations */
extern pop *g;
extern my_double *tt,*ttold;
#endif

#ifdef SALT
/* Populations */
extern pop *s;
extern my_double *ss,*ssold;
#endif

#ifdef TEMPERATURE_MELTING
extern my_double *ll,*llold, *hh;
#endif

#ifdef FLUID_RHEOLOGY
extern my_double *tau;
#endif














