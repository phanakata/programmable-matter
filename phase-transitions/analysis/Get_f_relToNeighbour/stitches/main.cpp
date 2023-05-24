#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "stdint.h"
#include "gsd.h"
#include "variables.h"
#include "gsd_read.h"
#include "analyze.h"
#include <string>

int NX,NY,RUN,afm,spx,spy, START;
double KAPPA;
int STEPS,FRAMES,CLONE;

int main(int argc, char **argv)
{
  FILE *therm,*lat;
  char filename[10];
  char init_strip[256],trajectory_file[256],thermalpos_file[256],lattice_file[256];
  int frame_cnt=0;
  double average_spins=0, average_spins2=0, average_spins4=0;
  
  /*   Output File     */
  FILE *bc;

   /*      Character array for directory pathname and filename     */
  char filepath[256];

  
   switch (argc){
     case 9:
       sscanf(argv[1],"%d",&NX);    
       sscanf(argv[2],"%d",&NY);
       sscanf(argv[3],"%s",filename);
       sscanf(argv[4],"%d",&START);
       sscanf(argv[5],"%d",&STEPS);
       sscanf(argv[6],"%d",&spx);
       sscanf(argv[7],"%d",&spy);
       sscanf(argv[8],"%d",&afm); 
       break;
     default:
       print_and_exit("Usage: %s NX NY filename START STEPS spx spy afm\n",
           argv[0]);
   }
   
   //START=0 get configurations from t=0
   //START=1 get configurations from t=1/4 of total frames
   //START=2 get configurations from t=1/2 of total frames
   
   FRAMES=STEPS/PERIOD;
   
   int frameStart = int(START*FRAMES/4);
   sprintf(trajectory_file, filename);
   frame_cnt=0;
   for(int frames=frameStart;frames<FRAMES;frames++)
     {
       load_gsd(trajectory_file,frames);
       average_spins=avg_spins();
       //slider_thermal+=spins;
       average_spins2 += average_spins*average_spins;
       average_spins4 += average_spins*average_spins*average_spins*average_spins;
       printf("%d\t%lf\n", frames, average_spins);
       frame_cnt++;
     }
   
 

   
   
   sprintf(filepath,"./binderCumulant%d", afm);
   //printf("Filename of Lattice Details: %s\n",filepath);
   bc = fopen(filepath, "w");
   //if (bc == NULL)
   //{
   //  print_and_exit("Could Not Open File:binderCumulant%d", afm);
   //}
   
   fprintf(bc,"%lf\t%lf", average_spins2/frame_cnt, average_spins4/frame_cnt);
   
   fclose(bc);
   

   return 0;
}
