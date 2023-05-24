#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "analyze.h"

/*	Slider average position	*/
double avg_Rg2(int frames)
{
   double slider_pos=0, min=1000, max=-1000, xx=0;
   double xcm=0, ycm=0, zcm=0, Rg2=0, Rg2_x=0, Rg2_y=0, Rg2_z=0;
   int slider_node=0;
   for(int i=0;i<N;i++)
   {   	
	xcm+=position[3*i];
	ycm+=position[3*i+1];
	zcm+=position[3*i+2];
	slider_node++;
	//int kk=1;
        //printf("%d\t%lf\t%lf\t%lf\n", i, position[i*3], position[i*3+1],position[i*3+2]);
        //if (position[3*i+kk]>max)
	//{
	//max=position[3*i+kk];
	//}
        //if (position[3*i+kk]<min) 
	//{
	//min = position[3*i+kk];
	//}

   }   
   xcm=xcm/slider_node;
   ycm=ycm/slider_node;
   zcm=ycm/slider_node;
   //printf("%d\t%lf\t%lf\t%lf\n", slider_node, xcm, ycm, zcm);
   for(int i=0; i<N; i++)
   {
	Rg2_x+=(position[3*i]-xcm)*(position[3*i]-xcm);
	Rg2_y+=(position[3*i+1]-ycm)*(position[3*i+1]-ycm);
        Rg2_z+=(position[3*i+2]-zcm)*(position[3*i+2]-zcm);

//+(position[3*i+1]-ycm)*(position[3*i+1]-ycm)+(position[3*i+2]-zcm)*(position[3*i+2]-zcm);  	 
       //Rg2+=(position[3*i+2]-zcm)*(position[3*i+2]-zcm);
   }	  
   Rg2_x = Rg2_x/slider_node;
   Rg2_y = Rg2_y/slider_node;
   Rg2_z = Rg2_z/slider_node;
   Rg2 = Rg2_x+ Rg2_y+Rg2_z;
		
   printf("%d\t%lf\t%lf\t%lf\t%lf\n", frames, Rg2, Rg2_x, Rg2_y, Rg2_z);
   return (Rg2);
}

