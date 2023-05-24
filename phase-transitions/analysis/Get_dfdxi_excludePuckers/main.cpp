#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "stdint.h"
#include "gsd.h"
#include "variables.h"
#include "gsd_read.h"
#include "analyze.h"

int NX, NY, RUN, CALC, DIM;
double DELTA;
int STEPS,FRAMES,CLONE;

int main(int argc, char **argv)
{
  FILE *therm,*lat, *fq2data;
  char init_strip[256],trajectory_file[256],thermalpos_file[256],lattice_file[256], fq2data_file[256];
  int frame_cnt=0;
  double h2_ave=0, h2_frame=0;
  
  switch (argc){
  case 6:
    sscanf(argv[1],"%d",&NX);    
    sscanf(argv[2],"%d",&NY);
    sscanf(argv[3],"%lf",&DELTA);
    sscanf(argv[4],"%d",&CALC);
    sscanf(argv[5],"%d",&STEPS);
    break;
     default:
       print_and_exit("Usage: %s NX NY DELTA CALC STEPS\n",
		      argv[0]);
  }

  
  FRAMES=STEPS/PERIOD;
  
  
  sprintf(trajectory_file,"traj.gsd");
  frame_cnt=0;
  double sum=0;
  double calc=0;
  for(int frames=FRAMES/2;frames<FRAMES;frames++)
    {
      load_gsd(trajectory_file,frames);
      //in each frame initialize and calculate sum (dfdx)^2 + (dfdy)^2
      //initialize_fmesh();
      //build_fmesh();
      if (CALC==1 or CALC==2)
	{
	  calc = build_fmesh_cd();
	  printf("%lf\n", calc);
	  sum += calc; //build_fmesh();//calculate_sum_dfdxi2();
	}
      else if (CALC==3 or CALC==4)
	{
	  calc = build_fmesh_cd_nodes(); //using nodes instead of square grid 
	  printf("%lf\n", calc);
	  sum += calc; //build_fmesh();//calculate_sum_dfdxi2();
	}
      else{
		printf("%lf\t%lf\n", find_xmax()-find_xmin(), find_ymax()-find_ymin());
	}
      frame_cnt++;
    }
  
  //printf("%lf", sum/frame_cnt);
  
  
  
  return 0;
}
