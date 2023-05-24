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
  int STEPS,TOTFRAMES,CLONE;

int main(int argc, char **argv)
{
  FILE *therm,*lat, *xyzdata, *spindata;
  char filename[10];
  char xyzdata_file[256], spin_file[256], init_strip[256],trajectory_file[256],thermalpos_file[256],lattice_file[256];
  int START;
  //double average_spins;
  
  switch (argc){
  case 9:
    sscanf(argv[1],"%d",&NX);    
    sscanf(argv[2],"%d",&NY);
    sscanf(argv[3],"%s",filename);
    sscanf(argv[4],"%d",&START);
    sscanf(argv[5],"%d",&TOTFRAMES);
    sscanf(argv[6],"%d",&spx);
    sscanf(argv[7],"%d",&spy);
    sscanf(argv[8],"%d",&afm); 
       break;
   default:
     print_and_exit("Usage: %s NX NY START FINAL spx spy afm\n",
		    argv[0]);
  }
 if (afm==1){
 	sprintf(xyzdata_file,"spins-staggered.xyz", afm);}
 else{
   sprintf(xyzdata_file,"spins.xyz", afm);}

 if (afm==1){
        sprintf(spin_file,"spinsConf-staggered.dat", afm);}
 else{
   sprintf(spin_file,"spinsConf.dat", afm);}


 xyzdata = fopen(xyzdata_file, "w");
  if (xyzdata == NULL)
    {
      print_and_exit("Could Not Open File to write xyz data");
    }
 
 
 spindata = fopen(spin_file, "w");
  if (spindata == NULL)
    {
      print_and_exit("Could Not Open File to write xyz data");
    }

initialize_xyz();     

sprintf(trajectory_file,filename);
int count_index=0;

  for (int frame=START; frame<TOTFRAMES; frame++)
{
  load_gsd(trajectory_file,frame);

  count_index = avg_spins(count_index);
 }
int ii =0; 
//planar PBC
int Nspins = (NX/2)*(NY/2);
//planar 
//int Nspins = (NX/2-2)*(NY/2-2);
//cylinder 
//int Nspins = (NX/2-1)*(NY/2);

int NspinsPos = 4*Nspins;
for (int frame=0; frame<TOTFRAMES-START; frame++)
{
 fprintf(xyzdata, "%d\n\n", Nspins); 
 for (int i=0; i<NspinsPos; i+=4){
 ii = i + NspinsPos*frame;
 //printf("%d\n", ii);
 fprintf (xyzdata, "%d %lf %lf %lf\n", int(xyz[ii]), xyz[ii+1]/2, xyz[ii+2]/2, xyz[ii+3]/2);
 fprintf (spindata, "%d\t", int(xyz[ii]));
}
 fprintf(spindata, "\n");
}
fclose(xyzdata);
fclose(spindata);

  return 0;
}
