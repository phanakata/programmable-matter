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
int TOTSTEPS,FRAMES,CLONE, START, NSAMPLE, DELTASTEPS;

int main(int argc, char **argv)
{
  FILE *therm,*lat, *fq2data;
  char init_strip[256],trajectory_file[256],thermalpos_file[256],lattice_file[256], fq2data_file[256];
  int frame_cnt=0;
  double h2_ave=0, h2_frame=0;
  
  switch (argc){
  case 8:
    sscanf(argv[1],"%d",&NX);    
    sscanf(argv[2],"%d",&NY);
    sscanf(argv[3],"%d",&START);
    sscanf(argv[4],"%d",&TOTSTEPS);
    sscanf(argv[5],"%d",&NSAMPLE);
    sscanf(argv[6],"%d",&SQUARED);
    sscanf(argv[7],"%d",&DIM);
    break;
     default:
       print_and_exit("Usage: %s NX NY START TOTSTEPS SQUARED DIM\n",
		      argv[0]);
  }
  // Avg, Height Squared ribbon profile path
  if (DIM==2)
	{
  	sprintf(fq2data_file,"fq2.dat");
  	}
  else
	{
	sprintf(fq2data_file,"fq2-1D.dat");
	}//printf("fq2 File: %s\n",fq2data_file);

  
  fq2data = fopen(fq2data_file, "w");
  if (fq2data == NULL)
    {
      print_and_exit("Could Not Open File to write fq2 data");
    }

  // Avg, Height Squared ribbon profile path
  //sprintf(fq2data_file,"fq2.dat");
  //printf("fq2 File: %s\n",fq2data_file);
  
  //START=START/PERIOD;//convert to frames 
  //TOTSTEPS = TOTSTEPS/PERIOD;
  DELTASTEPS=(TOTSTEPS-START)/NSAMPLE;
  
  sprintf(trajectory_file,"traj.gsd");
  frame_cnt=0;
  initialize_fq2(); // Initializing the fq2 array 
  for(int frames=START;frames<TOTSTEPS;frames+=DELTASTEPS)
    {
      
      load_gsd(trajectory_file,frames);
      if (SQUARED==0)
	{
	  h2_frame=avg_hgt();
	  printf("%lf\n", h2_frame);
	}
      else
	{
	  h2_frame=avg_hgt_sq();
	}
      h2_ave+=h2_frame;
     if (DIM==2)
	{ 
		calculate_fq2();
	}
      else
	{	
	      	calculate_fq2_1D();
	}
      //printf("%d\t%lf\n", frames, hcm);
      frame_cnt++;
    }
  //Average Height Fluctuation 
  if (DIM==2)
	{
	//int volume = NX*NX*NY*NY;
        int volume = (NX+1)*(NX+1)*(NY+1)*(NY+1);//OLD one 

  	avg_fq2(fq2data,frame_cnt);
  	double h2_ave_fromFourier = sum_fq2(); 
  	double sq_dfdx = sum_q2fq2();
  	fclose(fq2data);
  	printf("%lf\t%lf\t%lf\n", h2_ave/frame_cnt, h2_ave_fromFourier/(volume)/frame_cnt, sq_dfdx/(volume)/frame_cnt); 
	}
  else
	{
	int volume = NX*NX;
	avg_fq2_1D(fq2data,frame_cnt);
        double h2_ave_fromFourier = sum_fq2_1D();
        double sq_dfdx = sum_q2fq2_1D();
        fclose(fq2data);
        printf("%lf\t%lf\t%lf\n", h2_ave/frame_cnt, h2_ave_fromFourier/(volume)/frame_cnt, sq_dfdx/(volume)/frame_cnt);

	}
  //h2_ave/=frame_cnt;

  //initialize_fq2(); // Initializing the fq2 array 
  
  //for (int count_q=0; count_q<NX*NY; count_q++)
  //   {
  //    fq2[count_q]/=frame_cnt;
  //  }
  //printf("%lf\n", h2_ave);

  //FILE *fp;
  //fp = fopen("fq2.txt","w");
  //char lang[5][20] = {"C","C++","Java","Python","PHP"};
  //double qx=0, qy=0;
  //double delta_q = 2.*M_PI/(double)NX;

  
  //fprintf(fp,"Top 5 programming language\n");
  //int cont_q=0;
  // for (int nx=-NX/2; nx<NX/2; nx++)
  //  {
  //    for (int ny=-NY/2; ny<NY/2; ny++)
  //	{
  // qx=nx*delta_q;
  //	  qy=ny*delta_q;

	  //for (int i=0; i<5; i++)
  //	  fprintf(fp, "%lf\t%lf\t%lf\n", qx, qy, fq2[count_q]);
  //	}
  //}
  //fclose(fp);
  //return 0;
  
  
  
  return 0;
}
