#include "common_object.h"

void nusselt(int tstep, int flag){
  FILE *fout;
  char fname[128];
  int x, y, yp,ym;
  my_double Nusselt, nu, kappa;  
  my_double tmp1,tmp2,tmp3, tmp4, tmp5, vt, v2, t1, t2 , dyt;
  my_double vx_y[NY] ,  vy_y[NY];
  my_double vx2_y[NY], vy2_y[NY];
  my_double rho_y[NY];
  my_double nusselt_y[NY], vyt_y[NY] , dyt_y[NY];
  my_double t_y[NY], t2_y[NY];
  my_double  s_y[NY] , s2_y[NY];
  my_double lf_y[NY];

  nu =property.nu;
  kappa = property.kappa_t;


  if(flag==0){
    // fprintf(stderr, "tstep %d, flag %d...", tstep, flag);
    rbout.vx = rbout.vy = rbout.t = 0.0;
    rbout.vx2 = rbout.vy2 = rbout.t2 = 0.0;
    rbout.nusselt = rbout.dyt =  rbout.vyt = rbout.rho = 0.0;
#ifdef SALT
    rbout.s = rbout.s2 = 0.0;
#endif
#ifdef TEMPERATURE_MELTING
    rbout.lf = 0.0;
#endif

    for (y=0; y<NY; y++){
      nusselt_y[y] = 0.0;
      t_y[y] = vx_y[y] = vy_y[y] = 0.0;
      t2_y[y] = vx2_y[y] = vy2_y[y] = 0.0;
      rho_y[y] = 0.0;
#ifdef SALT
     s_y[y] = s2_y[y] = 0.0;
#endif
#ifdef TEMPERATURE_MELTING
      lf_y[y] = 0.0;
#endif

    }
  }

  for (y=1; y<NY+1; y++){
    yp=y+1;
    ym=y-1;
    for (x=1; x<NX+1; x++){
      
#ifdef METHOD_FORCING_GUO      
      tmp1 = (vx(p[IDX(y,x)]) + 0.5*force[IDX(y,x)].x)/m(p[IDX(y,x)]); 
      tmp2 = (vy(p[IDX(y,x)]) + 0.5*force[IDX(y,x)].y)/m(p[IDX(y,x)]); 
      //tmp1 = v[IDX(y,x)].vx;
      //tmp2 = v[IDX(y,x)].vy;
#else
      tmp1 = vx(p[IDX(y,x)])/m(p[IDX(y,x)]);
      tmp2 = vy(p[IDX(y,x)])/m(p[IDX(y,x)]);     
#endif

      tmp3 = t(g[IDX(y,x)]);
#ifdef SALT
      tmp4 = t(s[IDX(y,x)]);
#endif
#ifdef TEMPERATURE_MELTING
      tmp5 = ll[IDX(y,x)];
#endif
      rbout.vyt +=  tmp2*tmp3;  vyt_y[y-1] +=  tmp2*tmp3;
      rbout.vx2 +=  tmp1*tmp1;  vx2_y[y-1] +=  tmp1*tmp1;
      rbout.vy2 +=  tmp2*tmp2;  vy2_y[y-1] +=  tmp2*tmp2;
      rbout.t2  +=  tmp3*tmp3;  t2_y[y-1]  +=  tmp3*tmp3;
      rbout.vx +=  tmp1;        vx_y[y-1] +=  tmp1;
      rbout.vy +=  tmp2;        vy_y[y-1] +=  tmp2;
      rbout.t  +=  tmp3;        t_y[y-1]  +=  tmp3; 
      rbout.rho += m(p[IDX(y,x)]) ;        rho_y[y-1]  +=  m(p[IDX(y,x)]);   
      if(y==1){  rbout.dyt += (t(g[IDX(yp,x)]) - property.deltaT/2.0)/(1.5);  dyt_y[y-1] += (t(g[IDX(yp,x)]) - property.deltaT/2.0)/(1.5); }
      if(y==NY){ rbout.dyt += (-property.deltaT/2.0  - t(g[IDX(ym,x)]))/(1.5); dyt_y[y-1] += (-property.deltaT/2.0  - t(g[IDX(ym,x)]))/(1.5); }
      if(y>1 && y<NY){ rbout.dyt += 0.5*(t(g[IDX(yp,x)])-t(g[IDX(ym,x)])); dyt_y[y-1] += 0.5*(t(g[IDX(yp,x)])-t(g[IDX(ym,x)])); }

#ifdef SALT
      rbout.s  += tmp4;    s_y[y-1]  += tmp4;
      rbout.s2 += tmp4;    s2_y[y-1]  += tmp4*tmp4;
#endif
#ifdef TEMPERATURE_MELTING
      rbout.lf += tmp5; lf_y[y-1] += tmp5;
#endif
    }       

    //Nusselt += 1.0 + (vt/(my_double)NX)/(kappa*deltaT/(my_double)NY);  /* approximate computation of Nusselt number */

    //Nusselt_y[y-1]= (rbout.vyt/(my_double)NX - kappa*rbout.dyt/(my_double)NX)/(kappa*deltaT/(my_double)(NY)); 

    //Nusselt_y[y-1] = 1.0 + (vt/(my_double)NX)/(kappa*deltaT/(my_double)NY);  
  }

  if(flag==1){
    // fprintf(stderr, "flag %d, tstep %d\n", flag, tstep);
    rbout.nusselt = (rbout.vyt - property.kappa_t*rbout.dyt)/(property.kappa_t*property.deltaT/(my_double)NY);
    rbout.nusselt /= (2.0*NX*NY);
    rbout.vx2 /= 2.0*NX*NY;
    rbout.vy2 /= 2.0*NX*NY;
    rbout.t2 /= 2.0*NX*NY;
    rbout.t /= 2.0*NX*NY;
    rbout.rho /= 2.0*NX*NY;
    rbout.vx /= 2.0*NX*NY;
    rbout.vy /= 2.0*NX*NY;
#ifdef SALT
    rbout.s  /= 2.0*NX*NY;
    rbout.s2 /= 2.0*NX*NY;
#endif
#ifdef TEMPERATURE_MELTING
    rbout.lf /= 2.0*NX*NY; 
#endif
    for (y=1; y<NY+1; y++){
      dyt_y[y-1] /= 2.0*NX;
      vyt_y[y-1] /= 2.0*NX;
      vx_y[y-1] /= 2.0*NX;
      vy_y[y-1] /= 2.0*NX;
      t_y[y-1] /= 2.0*NX;
      nusselt_y[y-1] = (vyt_y[y-1] - property.kappa_t*dyt_y[y-1])/(property.kappa_t*property.deltaT/(my_double)NY); 
      rho_y[y-1] /= 2.0*NX;
      vx2_y[y-1] /= 2.0*NX;
      vy2_y[y-1] /= 2.0*NX;
      t2_y[y-1] /= 2.0*NX;
#ifdef SALT
      s_y[y-1]  /= 2.0*NX;
      s2_y[y-1] /= 2.0*NX;
#endif
#ifdef TEMPERATURE_MELTING
      lf_y[y-1] /= 2.0*NX; 
#endif
    }

    sprintf(fname,"nusselt.dat");
    fout = fopen(fname,"a");
    fprintf(fout,"%d %e %e %e %e %e %e %e %e\n",tstep, (double)rbout.nusselt, (double)rbout.vx2 , (double)rbout.vy2, (double)rbout.t2 , 
	    (double)rbout.vx, (double)rbout.vy, (double)rbout.t ,(double)rbout.rho);
    fclose(fout);

    sprintf(fname,"nusselt_y.dat");
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) fprintf(fout,"%d %e %e %e %e %e %e %e %e\n",y, (double)nusselt_y[y-1], (double)vx_y[y-1], (double)vy_y[y-1], (double)t_y[y-1], (double)rho_y[y-1],  (double)vx2_y[y-1], (double)vy2_y[y-1], (double)t2_y[y-1]);
    fclose(fout);
 

#ifdef SALT
    sprintf(fname,"salt.dat");
    fout = fopen(fname,"a");
    fprintf(fout,"%d %e %e\n",tstep, (double)rbout.s , (double)rbout.s2);
    fclose(fout);

    sprintf(fname,"salt_y.dat");
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) fprintf(fout,"%d %e %e\n",y, (double)s_y[y-1], (double)s2_y[y-1]);
    fclose(fout);
  
#endif

#ifdef TEMPERATURE_MELTING
    sprintf(fname,"melt.dat");
    fout = fopen(fname,"a");
    fprintf(fout,"%d %e\n",tstep, (double)rbout.lf);
    fclose(fout);

    sprintf(fname,"melt_y.dat");
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) fprintf(fout,"%d %e\n",y, (double)lf_y[y-1]);
    fclose(fout);
  
#endif
  }

}


/* print */
void print_fields(int tstep) 
{
  FILE *fout;
  char fname[128];
  int x, y, pp;

  if( tstep%time_dump_field == 0){

    /* Here dumps the populations */
#ifdef DUMP_POP
    for (pp=0; pp<9; pp++) {
      sprintf(fname,"%s/pop.%d.%d",OutDir,tstep,pp);
      fout = fopen(fname,"w");
      for (y=1; y<NY+1; y++) {
	for (x=1; x<NX+1; x++) 
	  fprintf(fout,"%d %d %g\n", x, y,
		  (double)p[IDX(y,x)].p[pp] );
	fprintf(fout,"\n");
      }
      fclose(fout);
    }
#endif


#ifdef FLUID
    /* Here dumps the velocity field */
    sprintf(fname,"%s/vel.%d",OutDir,tstep);
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) {
      for (x=1; x<NX+1; x++) 
	/*
#ifdef METHOD_FORCING_GUO
	fprintf(fout,"%d %d %g %g %g\n", x, y,
		(v[IDX(y,x)].vx + 0.5*force[IDX(y,x)].x)/dens[IDX(y,x)],
		(v[IDX(y,x)].vy + 0.5*force[IDX(y,x)].y)/dens[IDX(y,x)], dens[IDX(y,x)]);  
#else
	*/
	fprintf(fout,"%d %d %g %g %g\n", x, y,
		(double)v[IDX(y,x)].vx,
		(double)v[IDX(y,x)].vy, (double)dens[IDX(y,x)]);  
      /*
#endif  
      */    
      fprintf(fout,"\n");
    }
    fclose(fout);
#endif

#ifdef TEMPERATURE
    sprintf(fname,"%s/temp.%d",OutDir,tstep);
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) {
      for (x=1; x<NX+1; x++) 
	fprintf(fout,"%d %d %g\n", x, y, (double)tt[IDX(y,x)] );
      fprintf(fout,"\n");
    }
    fclose(fout);
#endif

#ifdef TEMPERATURE_MELTING
    /* Here dumps the liquid fraction field */
    sprintf(fname,"%s/melt.%d",OutDir,tstep);
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) {
      for (x=1; x<NX+1; x++) 
	fprintf(fout,"%d %d %g %g %g\n", x, y, (double)ll[IDX(y,x)],(double)(ll[IDX(y,x)]-llold[IDX(y,x)]), (double)hh[IDX(y,x)]);
      fprintf(fout,"\n");
    }
    fclose(fout);
#endif

#ifdef SALT
    /* Here dumps the liquid fraction field */
    sprintf(fname,"%s/salt.%d",OutDir,tstep);
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) {
      for (x=1; x<NX+1; x++) 
	fprintf(fout,"%d %d %g\n", x, y, (double)ss[IDX(y,x)]);
      fprintf(fout,"\n");
    }
    fclose(fout);
#endif

  }/*if tstep*/

}


/* print */
void print_fields_h5(int tstep) 
{
  FILE *fout;
  char fname[128];
  int x, y, pp;

  if( tstep%time_dump_field == 0){

    /* Here dumps the populations */
#ifdef DUMP_POP
    for (pp=0; pp<9; pp++) {
      sprintf(fname,"%s/pop.%d.%d",OutDir,tstep,pp);
      fout = fopen(fname,"w");
      for (y=1; y<NY+1; y++) {
	for (x=1; x<NX+1; x++) 
	  fprintf(fout,"%d %d %g\n", x, y,
		  (double)p[IDX(y,x)].p[pp] );
	fprintf(fout,"\n");
      }
      fclose(fout);
    }
#endif


#ifdef FLUID
    /* Here dumps the velocity field */
    sprintf(fname,"%s/vel.%d",OutDir,tstep);
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) {
      for (x=1; x<NX+1; x++) 
	fprintf(fout,"%d %d %g %g %g\n", x, y,
	  (double)(v[IDX(y,x)].vx/dens[IDX(y,x)]),
	  (double)(v[IDX(y,x)].vx/dens[IDX(y,x)]), (double)dens[IDX(y,x)]);        
      fprintf(fout,"\n");
    }
    fclose(fout);
#endif

#ifdef TEMPERATURE
    sprintf(fname,"%s/temp.%d",OutDir,tstep);
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) {
      for (x=1; x<NX+1; x++) 
	fprintf(fout,"%d %d %g\n", x, y, (double)tt[IDX(y,x)] );
      fprintf(fout,"\n");
    }
    fclose(fout);
#endif

#ifdef TEMPERATURE_MELTING
    /* Here dumps the liquid fraction field */
    sprintf(fname,"%s/melt.%d",OutDir,tstep);
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) {
      for (x=1; x<NX+1; x++) 
	fprintf(fout,"%d %d %g %g %g\n", x, y, (double)ll[IDX(y,x)],(double)ll[IDX(y,x)]-(double)llold[IDX(y,x)], (double)hh[IDX(y,x)]);
      fprintf(fout,"\n");
    }
    fclose(fout);
#endif

#ifdef SALT
    /* Here dumps the liquid fraction field */
    sprintf(fname,"%s/salt.%d",OutDir,tstep);
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) {
      for (x=1; x<NX+1; x++) 
	fprintf(fout,"%d %d %g\n", x, y, (double)ss[IDX(y,x)]);
      fprintf(fout,"\n");
    }
    fclose(fout);
#endif

  }/*if tstep*/

}
