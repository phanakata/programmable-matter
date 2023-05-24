#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "stdint.h"
#include "gsd.h"
#include "variables.h"
#include "gsd_read.h"
#include "analyze.h"

  int NX,NY,RUN,afm,spx,spy, vis;
  double KAPPA;
  int STEPS,TOTFRAMES,CLONE;

int main(int argc, char **argv)
{
  FILE *therm,*lat, *xyzdata;
  char filename[10];
  char xyzdata_file[256], init_strip[256],trajectory_file[256],thermalpos_file[256],lattice_file[256];
  int START;
  //double average_spins;
  
  switch (argc){
  case 10:
    sscanf(argv[1],"%d",&NX);    
    sscanf(argv[2],"%d",&NY);
    sscanf(argv[3],"%s",filename);
    sscanf(argv[4],"%d",&START);
    sscanf(argv[5],"%d",&TOTFRAMES);
    sscanf(argv[6],"%d",&spx);
    sscanf(argv[7],"%d",&spy);
    sscanf(argv[8],"%d",&afm);
    sscanf(argv[9],"%d",&vis);  
       break;
   default:
     print_and_exit("Usage: %s NX NY filename STARTFRAME FINALFRAME spx spy afm visualizespins\n",
		    argv[0]);
  }
 if (vis==1){
 if (afm==1){
 	sprintf(xyzdata_file,"spinsPBC-staggered.xyz", afm);}
 else{
   sprintf(xyzdata_file,"spinsPBC.xyz", afm);}


 xyzdata = fopen(xyzdata_file, "w");
  if (xyzdata == NULL)
    {
      print_and_exit("Could Not Open File to write xyz data");
    }
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
int Nspins = (NX/spx)*(NY/spy);

//pzh: to be deleted 
//planar 
//int Nspins = (NX/2-2)*(NY/2-2);
//cylinder 
//int Nspins = (NX/2-1)*(NY/2);

int NspinsPos = 4*Nspins;
double morder=0.;
int i_right=0, i_up=0, NXSpins=NX/spx;
for (int frame=0; frame<TOTFRAMES-START; frame++)
{ 
 morder=0.;
 if (vis==1){fprintf(xyzdata, "%d\n\n", Nspins);}
 for (int i=0; i<NspinsPos; i+=4){
 ii = i + NspinsPos*frame;
 //printf("%d\n", ii);
 if (vis==1)
	{fprintf (xyzdata, "%d %lf %lf %lf\n", int(xyz[ii]), xyz[ii+1]/2, xyz[ii+2]/2, xyz[ii+3]/2);}
 
 //right side 
 i_right = ii + 4;
 //wrap x-pbc 
 if (i_right%(4*NXSpins)==0)
	{
	i_right=i_right-4*NXSpins; 	
	} 

 morder += ((xyz[ii]-0.5)*2)*((xyz[i_right]-0.5)*2); //need to covnert +1/-1
 //up side
 i_up = ii+4*NXSpins;
 //wrap y-pbc
 if (i_up>=NspinsPos*(frame+1))
	{
        i_up=i_up-NspinsPos;
	}
 morder +=((xyz[ii]-0.5)*2)*((xyz[i_up]-0.5)*2);
 //debuging 
 //printf("%d\t%d\t%d\t%lf\t%lf\n", ii, i_right, i_up, ((xyz[ii]-0.5)*2)*((xyz[i_right]-0.5)*2), ((xyz[ii]-0.5)*2)*((xyz[i_up]-0.5)*2));
 }
 printf("%lf\n", morder/(double)Nspins/2.);
}
if (vis==1){
fclose(xyzdata);
}
  return 0;
}
