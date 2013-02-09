#include "common_object.h"




 /*********************************************                                                               
 *    y                                       *
 *    ^                                       *
 *    |                                       *
 *    |                                       *
 *    |                                       *
 *    +----------> x                          *
 *                                            *                                           
 *         x  y                               *
 *  0    (+0,+0)         6    2    5          *
 *  1    (+1,+0)          \   |   /           *
 *  2    (+0,+1)           \  |  /            *
 *  3    (-1,+0)            \ | /             *
 *  4    ( 0,-1)      3 <---- 0 ----> 1       *
 *  5    (+1,+1)            / | \             *
 *  6    (-1,+1)           /  |  \            *
 *  7    (-1,-1)          /   |   \           *
 *  8    (+1,-1)         7    4    8          * 
 *                                            *
 *********************************************/

/* attempt to  write free surface flow routine */

#ifdef FLUID_SURFACE
void surface(int i){
  int x, y , pp;
  int dx, dy;
  pop dm;

  /* store previous mass fraction */
  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++) {
      mmold[IDX(y,x)] = mm[IDX(y,x)]; 
    }

  /* compute mass fraction in each cell*/
  
  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++) 
       for (pp=0; pp<9; pp++){ 
	dx=cx[pp];
	dy=cy[pp];
	cell_type=mm[IDX(y+dy,x+dx)];

	if(cell_type==0){
	  dm[pp] = 0.0;
	} else if(cell_type==1){
	  dm[pp] = p[IDX(y+dy,x+dx)].[invp[pp]] - p[IDX(y,x)].[pp];
	} else {
	  dm[pp] = 0.5*(mm[IDX(y,x)]+mm[IDX(y+dy,x+dx)])*(p[IDX(y+dy,x+dx)].[invp[pp]] - p[IDX(y,x)].[pp]);
	}


	


	

      }


}
#endif
