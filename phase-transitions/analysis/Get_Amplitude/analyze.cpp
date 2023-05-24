#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "analyze.h"

/*	Slider average position	*/
double avg_h2()
{
  double hcm=0, h2=0, min=1000, max=-1000, xx=0;
  int tot_node=0;
  for(int i=0;i<N;i++)
    {          
      //xx=position[3*i];
      //pzh gsd reader doesn't recognize different particle ID , all particle ID=0 
      //if(particleID[i]==0)
      //if (i==N/2)
	//{
	//printf ("%lf\n", position[3*i+2]);
	//}      
      hcm+=position[3*i+2];
      tot_node++;
      
      //if (xx>max)
	//{
	 // max =xx;
	//}
      
      //if (xx<min)
	//{
	  //min =xx;
	//}
      
     
    }
  
  //printf ("%lf\n", max-min);
  hcm = hcm/tot_node;
  //printf
  //printf ("hcm: %lf\t%lf\n", hcm, hcm/tot_node);
 
  for(int i=0;i<N;i++)
    {          
      
      //xx=position[3*i];
      //pzh gsd reader doesn't recognize different particle ID , all particle ID=0 
      //if(particleID[i]==0)
      
      xx=(position[3*i+2]-hcm); //*(position[3*i+2]-hcm);
      if (xx>max)
	{
	max =xx;
	} 
      if (xx<min)
	{
	min =xx;
	}

//tot_node++;
      //hcm=hcm/(double)tot_node;
    }
   

   
   
   //printf("%d\t", slider_node);
   printf ("%lf\t%lf\n", max-min, hcm);
   return (max-min);
}


double avg_hgt_sq()
{
   double slider_pos=0, min=1000, max=-1000, xx=0;
   int slider_node=0;
   for(int i=0;i<N;i++)
   {   	//printf("Hi");
	xx=position[3*i];
	//pzh gsd reader doesn't recognize different particle ID , all particle ID=0 
	if(particleID[i]==0)
	{	
		slider_pos+=position[3*i+2]*position[3*i+2];
		slider_node++;
	
		if (xx>max)
		{
			max =xx;
		}
	
		if (xx<min)
		{
			min =xx;
		}
	
	}
   }
   //printf("%d\t", slider_node);
   printf ("%lf\t", max-min);
   return (slider_pos/slider_node);
}

