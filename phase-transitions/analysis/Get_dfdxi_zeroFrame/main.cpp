#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "stdint.h"
#include "gsd.h"
#include "variables.h"
#include "gsd_read.h"
#include "analyze.h"

int NX, NY, RUN, SQUARED, DIM;
double KAPPA;
int STEPS,FRAMES,CLONE;

int main(int argc, char **argv)
{
  FILE *therm,*lat, *fq2data;
  char init_strip[256],trajectory_file[256],thermalpos_file[256],lattice_file[256], fq2data_file[256];
  int frame_cnt=0;
  double h2_ave=0, h2_frame=0;
  
  switch (argc){
  case 5:
    sscanf(argv[1],"%d",&NX);    
    sscanf(argv[2],"%d",&NY);
    sscanf(argv[3],"%lf",&KAPPA);
    sscanf(argv[4],"%d",&STEPS);
    break;
     default:
       print_and_exit("Usage: %s NX NY KAPPA STEPS\n",
		      argv[0]);
  }

  
  FRAMES=STEPS/PERIOD;
  
  
  sprintf(trajectory_file,"traj.gsd");
  frame_cnt=0;
  double sum=0;
  for(int frames=FRAMES/2;frames<FRAMES;frames++)
    {
      load_gsd(trajectory_file,frames);
      //in each frame initialize and calculate sum (dfdx)^2 + (dfdy)^2
      initialize_fmesh();
      build_fmesh();
      sum +=calculate_sum_dfdxi2();
      frame_cnt++;
    }
  
  printf("%lf", sum/frame_cnt);
  
  
  
  return 0;
}
