#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "stdint.h"
#include "gsd.h"
#include "variables.h"
#include "gsd_read.h"
#include "analyze.h"

  int NX,NY,RUN,afm,spx,spy;
  double KAPPA;
  int STEPS,FRAMES,CLONE;

int main(int argc, char **argv)
{
  FILE *therm,*lat;
  char init_strip[256],trajectory_file[256],thermalpos_file[256],lattice_file[256];
  int FRAME;
  double average_spins;
  
  switch (argc){
  case 9:
    sscanf(argv[1],"%d",&NX);    
    sscanf(argv[2],"%d",&NY);
    sscanf(argv[3],"%lf",&KAPPA);
    sscanf(argv[4],"%d",&FRAME);
    sscanf(argv[5],"%d",&STEPS);
    sscanf(argv[6],"%d",&spx);
    sscanf(argv[7],"%d",&spy);
    sscanf(argv[8],"%d",&afm); 
       break;
   default:
     print_and_exit("Usage: %s NX NY KAPPA FRAME STEPS spx spy afm\n",
		    argv[0]);
  }
  sprintf(trajectory_file,"traj.gsd");
  //get spin configurations for a specific frame
  load_gsd(trajectory_file,FRAME);
  average_spins=avg_spins();
  
  return 0;
}
