#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "stdint.h"
#include "gsd.h"
#include "variables.h"
#include "gsd_read.h"
#include "analyze.h"

int NX,NY,RUN, SQUARED;
double KAPPA;
int STEPS,FRAMES,CLONE;

int main(int argc, char **argv)
{
  FILE *therm,*lat;
  char init_strip[256],trajectory_file[256],thermalpos_file[256],lattice_file[256];
  int frame_cnt=0;
  double h2_ave=0, h2_frame=0;
  
  switch (argc){
  case 6:
    sscanf(argv[1],"%d",&NX);    
    sscanf(argv[2],"%d",&NY);
    sscanf(argv[3],"%lf",&KAPPA);
    sscanf(argv[4],"%d",&SQUARED);
    sscanf(argv[5],"%d",&STEPS);
    //sscanf(argv[6],"%d",&CLONE); 
    break;
     default:
       print_and_exit("Usage: %s NX NY KAPPA SQUARED STEPS\n",
		      argv[0]);
  }
  
  FRAMES=STEPS/PERIOD;
  
  
  sprintf(trajectory_file,"traj.gsd");
  frame_cnt=0;
  for(int frames=FRAMES/2;frames<FRAMES;frames++)
    {
      load_gsd(trajectory_file,frames);
      if (SQUARED==0)
	{
	  h2_frame=avg_h2();
	  printf("%lf\n", h2_frame);
	}
      else
	{
	  h2_frame=avg_hgt_sq();
	}
      h2_ave+=h2_frame;
      //printf("%d\t%lf\n", frames, hcm);
      frame_cnt++;
    }
  h2_ave/=frame_cnt;
  //printf("%lf\n", h2_ave);
	 
  
  
  return 0;
}
