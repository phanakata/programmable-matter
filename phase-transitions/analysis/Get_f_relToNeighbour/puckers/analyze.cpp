#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "analyze.h"

/*	Average spin	*/
double avg_spins()
{
   double slider_pos=0, min=1000, max=-1000, xx=0;
   double total_spins=0;
   double spin; 
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
       for(int i=edgex+spx; i<NX-spx;i+=spx)
	 for (int j=edgey+spy; j<NY-spy; j+=spy)
	   {
	     kk++;
	     //staggered operator 
	     if (afm==1 or afm==3)
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
		 //B relative A (MIDPOINT)
		 rBA[jj]=(position[3*(j*NX+i+1)+jj]-position[3*(j*NX+i-1)+jj])/2.; //use mid point 
		 //C relative to A
		 rCA[jj]=position[3*(j*NX+i+NX)+jj]-position[3*(j*NX+i-1)+jj];
		 //defect position relative to A
		 rXA[jj]=position[3*(j*NX+i)+jj]-position[3*(j*NX+i-1)+jj];
		 
		 //relative pos X respect to the mid point 
		 delta_z[jj]=rXA[jj]-rBA[jj];
	       } 
	     
	     //square plane formed by 3 NNN

	    
	     
	     rCA2=0.;
	     rBA2=0.;
	     //rXA2=0.;
	     delta_z2=0.;
	     //get absolute value 
	     for (int jj=0; jj<3; jj++)
	       {
		 rBA2 +=rBA[jj]*rBA[jj]; 
		 rCA2 +=rCA[jj]*rCA[jj];
                 //rXA2 +=rXA[jj]*rXA[jj];
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
	   
	     
	     // we can estimate the 'height' relative to the plane as f = sqrt (rXA^2 -	 

	     
	     if (afm==1 or afm==0){			   
	     if(dot_product>0)
		       
	       {
		 spin=delta_z2;
		 total_spins+= double(spin*stg);
		 cnt_spins++;
	       }
	     
	     else
	       {
		 spin=-delta_z2;
		 total_spins+= double(spin*stg);
		 cnt_spins++;
	       }
		}
	     else if (afm==2){

		total_spins+=delta_z2; //absolute |m|
                cnt_spins++; 
			}

              else if (afm==3){

                total_spins+= stg*position[3*(j*NX+i)+2]; //relative to zero plane
                cnt_spins++;
                        }
	     else if (afm==4){

		total_spins+= position[3*(j*NX+i)+2]*position[3*(j*NX+i)+2]; //relative to zero plane
                cnt_spins++;


		}
		else{
		total_spins+=delta_z2*delta_z2; //absolute
                cnt_spins++;

		}


             //if (spin*stg!=1)
		//	{
	     //printf("%d ", spin);
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
   //printf("\nTOTAL SPINS: %d\n", cnt_spins);
   

   
   
   
  
   //printf ("%lf\t", max-min);
   return (total_spins/cnt_spins);
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
