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
  double slider_thermal=0, hcm=0;
  
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
	  hcm=avg_hgt();
	}
      else
	{
	  hcm=avg_hgt_sq();
	}
      slider_thermal+=hcm;
      printf("%d\t%lf\n", frames, hcm);
      frame_cnt++;
    }
  slider_thermal/=frame_cnt;
  
  
  
  return 0;
}
