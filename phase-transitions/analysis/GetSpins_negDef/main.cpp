/*
 * Paul Hanakata 2020-2022 - Harvard University 
 * C++ code to calculate average magnetization, staggered magentization and Binder cummulants
 * The "spin" was measured by digitizing the height of the buckled unit relative its local neighbours 
 * */

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

int NX,NY,afm,spx,spy;
double KAPPA;
int TOTSNAPSHOTS,START;

int main(int argc, char **argv)
{
  char filename[10];
  char trajectory_file[256];
  int frame_cnt=0;
  double average_spins=0, average_spins2=0, average_spins4=0;
  
  /*   Output File     */
  FILE *bc;

   /*      Character array for directory pathname and filename     */
  char filepath[256];

  
   switch (argc){
     case 9:
       sscanf(argv[1],"%d",&NX); //number of sites along x     
       sscanf(argv[2],"%d",&NY); //number of sites along y
       sscanf(argv[3],"%s",filename); //file name
       sscanf(argv[4],"%d",&START); //start (in units of frame/snapshot) 
       sscanf(argv[5],"%d",&TOTSNAPSHOTS); //end (in units of frame/snaphsot)
       sscanf(argv[6],"%d",&spx); //spacing between buckled unit in x 
       sscanf(argv[7],"%d",&spy); //spacing between buckled unit in y 
       sscanf(argv[8],"%d",&afm); //apply staggered (1) or not (0)
       break;
     default:
       print_and_exit("Usage: %s NX NY filename START TOTSNAPSHOTS spx spy afm\n",
           argv[0]);
   }
   
   sprintf(trajectory_file, filename);
   frame_cnt=0;
   for(int frames=START;frames<TOTSNAPSHOTS;frames++)
     {
       load_gsd(trajectory_file,frames);
       average_spins=avg_spins();
       average_spins2 += average_spins*average_spins;
       average_spins4 += average_spins*average_spins*average_spins*average_spins;
       printf("%d\t%lf\n", frames, average_spins);
       frame_cnt++;
     }
   
   sprintf(filepath,"./binderCumulant%d", afm);
   bc = fopen(filepath, "w");
   
   fprintf(bc,"%lf\t%lf", average_spins2/frame_cnt, average_spins4/frame_cnt);
   
   fclose(bc);
   

   return 0;
}
