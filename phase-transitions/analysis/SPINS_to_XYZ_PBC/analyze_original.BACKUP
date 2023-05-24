#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "analyze.h"

/*	Slider average position	*/
double xyz[NMAX*4*100];

int initialize_xyz()
{ 
  //planar
  int NN = (NX/2-2)*(NY/2-2)*4*100;
  //cylinder
  //int NN = (NX/2-1)*(NY/2)*100/4;	
  for(int i=0;i<NN; i++)
    {
      xyz[i]=0;
    }
  return 0;
}


int avg_spins(int count_index)
{
   //double slider_pos=0, min=1000, max=-1000, xx=0;
   double total_spins=0;
   int spin; 
   int cnt_spins=0;
   int pbc, edgex=spx-2, edgey=spy-2, stg;
   double rBA[3], rCA[3], rXA[3], delta_z[3]; //X=defect
   double rBA_cross_rCA[3];
   double rBA2=0, rCA2=0, delta_z2=0, dot_product=0, dz;
   int kk=0;	
   //we exclude boundaries for now
   pbc=1;//pbc true
   if (pbc==1) //periodic 
     {
//       for(int i=edgex+spx; i<NX;i+=spx)
//         for (int j=edgey; j<NY; j+=spy) 

//planr

    for(int i=edgex+spx; i<NX-spx;i+=spx)
	 for (int j=edgey+spy; j<NY-spy; j+=spy)
	   {
	     kk++;
	     //staggered operator 
	     if (afm==1)
	       {
		 stg = pow(-1, (i-edgex)/spx+(j-edgey)/spy);
	       }
	     else
	       {
		 stg=1;
	       }
	     //most simple triangle plane formed by 3 NN (top quadrant)
	     for (int jj=0; jj<3; jj++)
	       {
		 //B relative A
		 rBA[jj]=(position[3*(j*NX+i+1)+jj]-position[3*(j*NX+i-1)+jj])/2.; //use mid point 
		 //C relative to A
		 rCA[jj]=position[3*(j*NX+i+NX)+jj]-position[3*(j*NX+i-1)+jj];
		 //defect position relative to A
		 rXA[jj]=position[3*(j*NX+i)+jj]-position[3*(j*NX+i-1)+jj];
		 
		 //relative pos respect to the mid point 
		 delta_z[jj]=rXA[jj]-rBA[jj];
	       } 
	     
	     //square plane formed by 3 NNN

	    
	     
	     rCA2=0.;
	     rBA2=0.;
	     delta_z2=0.;
	     //get absolute value 
	     for (int jj=0; jj<3; jj++)
	       {
		 rBA2 +=rBA[jj]*rBA[jj]; 
		 rCA2 +=rCA[jj]*rCA[jj];
		 delta_z2 +=delta_z[jj]*delta_z[jj];
	       }
	     
	     rBA2 = sqrt(rBA2);
	     rCA2 = sqrt(rCA2);
	     delta_z2 = sqrt(delta_z2);
	     
	     //cross product 
	     //rBA cross rCA in unit vector! 
	     rBA_cross_rCA[0] = (rBA[1]*rCA[2]-rCA[1]*rBA[2])/rBA2/rCA2;
	     rBA_cross_rCA[1] = -(rBA[0]*rCA[2]-rCA[0]*rBA[2])/rBA2/rCA2;
	     rBA_cross_rCA[2] = (rBA[0]*rCA[1]-rCA[0]*rBA[1])/rBA2/rCA2;
	     
	     //printf("nvec: %lf %lf %lf\n", rBA_cross_rCA[0],  rBA_cross_rCA[1], rBA_cross_rCA[2]); 
	     //printf("nvec2: %lf %lf %lf\n", delta_z[0]/delta_z2,  delta_z[1]/(delta_z2), delta_z[2]/(delta_z2));      
	     //dot product 
	     dot_product=0;
	     for (int jj=0; jj<3; jj++)
	       {
		 dot_product+=rBA_cross_rCA[jj]*delta_z[jj]/(delta_z2);
	       }
	     
	     
	     //printf("\nID: %d\n", particleID[(j*NX+i)]);
	   
	     
	     

	     
		   
	     if(dot_product>0)
		       
	       {
		 spin=1;
                 xyz[count_index]=(spin*stg+1)/2; 
                 count_index++;
                 xyz[count_index]=position[3*(j*NX+i)];
                 count_index++;
                 xyz[count_index]=position[3*(j*NX+i)+1];
		 count_index++;
                 xyz[count_index]=position[3*(j*NX+i)+2];	
		 count_index++;
                 //total_spins+= double(spin*stg);
		 cnt_spins++;
	       }
	     
	     else
	       {
		 spin=-1;
		 xyz[count_index]=(spin*stg+1)/2; 
                 count_index++;
                 xyz[count_index]=position[3*(j*NX+i)];
                 count_index++;
                 xyz[count_index]=position[3*(j*NX+i)+1];
                 count_index++;
                 xyz[count_index]=position[3*(j*NX+i)+2];		 
		 count_index++;
                 spin=-1;
                 total_spins+= double(spin*stg);
		 cnt_spins++;
	       }
             //if (spin*stg!=1)
		//	{
	     //printf("%d ", stg*spin);
	     //}//total_spins+= double(spin);
	     //cnt_spins++;
	   }
     }
   else
     {
       
       //test 
       for(int i=edgex; i<NX;i+=spx)
	 for (int j=edgey; j<NY; j+=spy)
	   {
	     kk++;
	     //staggered operator 
	     if (afm==1)
	       {
		 stg = pow(-1, (i-edgex)/spx+(j-edgey)/spy);
	       }
	     else
	       {
		 stg=1;
	       }
	     
	     
	     
		
	     dz=position[3*(j*NX+i)+2];
	     //printf("\nID: %d\n", particleID[(j*NX+i)]);
	     //bumps
	     //
	     if (dz>position[3*(j*NX+i+NX)+2] and dz>position[3*(j*NX+i-NX)+2] and dz>position[3*(j*NX+i-1)+2] and dz>position[3*(j*NX+i+1)+2] ) 
	       
	       {
		 spin=1;
		 total_spins+= double(spin*stg);
		 cnt_spins++;
	       }
	     
	     else
	       {
		 spin=-1;
		 total_spins+= double(spin*stg);
		 cnt_spins++;
	       }
	      

	    

	   }
       //will add for nonPBC
     }
   //printf("\nTOTAL SPINS: %d\t%d\n", cnt_spins, kk);
   

   
   
   
  
   //printf ("%lf\t", max-min);
   return count_index;
}

/*

 int kk=0; 
  for(int i=3; i<NX-3;i+=spx)
    for (int j=3; j<NY-3; j+=spy)
      {
	kk++;
	dz=position[3*(j*NX+i)+2];
	//printf("\nID: %d\n", particleID[(j*NX+i)]);
	//bumps
	//
	if (dz>position[3*(j*NX+i+NX)+2] and dz>position[3*(j*NX+i-NX)+2] and dz>position[3*(j*NX+i-1)+2] and dz>position[3*(j*NX+i+1)+2] )        
	  {
	    spin=1;
	    total_spins+= double(spin);
	    cnt_spins++;
	  }
        
	else
	  {
	    spin=-1;
	    total_spins+= double(spin);
	    cnt_spins++;
	  }
	//printf("%d ", spin);
	//total_spins+= double(spin);
	//cnt_spins++;
      }
  //printf("\nTOTAL SPINS: %d\t%d\n", cnt_spins, kk);
  */
     //     if (dz>position[3*(j*NX+i+NX)+2] and dz>position[3*(j*NX+i-NX)+2] and dz>position[3*(j*NX+i-1)+2] and dz>position[3*(j*NX+i+1)+2] ) 
