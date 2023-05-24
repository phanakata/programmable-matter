#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "analyze.h"

//double fRq=0, fIq=0, qx=0, qy=0;
double fq2[NMAX]; 


double avg_hgt()
{
  double hgt=0;
  int node_cnt=0;
  for(int i=0;i<N;i++)
    {
      //if(particleID[i]==0 || particleID[i]==4)
      //{
      hgt+=position[3*i+2];
      node_cnt++;
      //}
    }
  return (hgt/node_cnt);
}

/*	Average <h^2> = 1/N * Sum_i(z_i - <z>)^2	*/
double avg_hgt_sq()
{
  double hgtSq=0;
  double h_avg = avg_hgt();
  int node_cnt=0;
  for(int i=0;i<N;i++)
    {
      //if(particleID[i]==0 || particleID[i]==4)
      //{
      hgtSq+=pow((position[3*i+2]-h_avg),2);
      node_cnt++;
      //}
    }
  return (hgtSq/node_cnt);
   
}

double get_LX()
{
  double LMAX=-100000, LMIN=10000, xx=0.;
  int node_cnt=0;
  for(int i=0;i<N;i++)
    {
      xx = position[i];
      if (xx>LMAX)
	{
	  LMAX=xx;
	}
      if (xx<LMIN)
	{
	  LMIN=xx;
	}
    }
  return (LMAX-LMIN);

}




int initialize_fq2()
{
   
  for(int i=0; i<(NX); i++)
    {
      fq2[i]=0;
    }
  return 0;
}


double sum_fq2()
{
  double delta_q = 2.*M_PI/(double)NX;
  double sum=0; 
  for(int i=0; i<NX/2; i++)
    {
      sum += (2.*M_PI)*(double)(1+i)*fq2[i]; 		
      //sum += (i+1)fq2[i]*(double)NY; //need to multipl by NY (#number of directios)		


    }
  
  return sum;
}

double sum_fq2_1D()
{
  double sum=0; 
  for(int i=0; i<NX; i++)
    {
      sum +=fq2[i];
    }
  
  return sum;
}


double sum_q2fq2()
{
  double sum=0, q=0;
  double delta_q = 2.*M_PI/(double)NX;
  //2.*M_PI/get_LX(); //2.*M_PI/(double)NX;
  int nx=0, ny=0, count_q=0;
  for (nx=0; nx<NX/2; nx++)
    {
      
      q=(double) (nx+1);
      
	  
      sum +=(2.*M_PI*q)*(q*q)*fq2[count_q];
      count_q++;
      
    }
  return sum;
}



double sum_q2fq2_1D()
{
  double sum=0, qx=0., qy=0.;
  double delta_q = 2.*M_PI/(double)NX;
  int nx=0, ny=0, count_q=0;
  //for (nx=-NX/2; nx<NX/2; nx++)
  for (nx=0; nx<NX/2; nx++)
    {
        qx=delta_q * (double) nx;
        sum +=(qx*qx)*fq2[count_q];
        count_q++;
    }
  return sum;
}


int calculate_fq2_1D()
{
  double delta_q = 2.*M_PI/(double)NX;
  double fRq=0, fIq=0, qx=0., qy=0., xx=0., yy=0., zz=0.;
  int count_q =0, nx=0, ny=0;
  double hcm=avg_hgt();
  for (nx=-NX/2; nx<NX/2; nx++)
    {
	  qx=delta_q * (double) nx;
	  fRq=0.;
          fIq=0.; 
	  for(int i=0;i<N;i++)
	    {          
	      xx = position[3*i];
	      zz = position[3*i+2]-hcm;
	      
	      fRq += zz * cos (qx * xx);
	      fIq += zz * sin (qx * xx);
	      
	    }
	  fq2[count_q]+=fRq*fRq + fIq*fIq;   
          count_q++;
	}
  return 0;
}
int calculate_fq2_new()
{
  double delta_q = 2.*M_PI/(double)NX;
  double fRq=0, fIq=0, qx=0., qy=0., xx=0., yy=0., zz=0.;
  int count_q =0, nx=0, ny=0;
  //double sum_fq2=0, fq2_abs=0;
  double hcm=avg_hgt();
  for (nx=0; nx<NX; nx++)
    {
      for (ny=0; ny<NY; ny++)
        {
          qx=delta_q * (double) nx;
          qy=delta_q * (double) ny;

          fRq=0.;
          fIq=0.;
          for(int i=0;i<N;i++)
            {
              xx = position[3*i]+double(NX)/2;
              yy = position[3*i+1]+double(NY)/2;
              zz = position[3*i+2]-hcm;


              fRq += zz * cos (qx * xx + qy * yy);
              fIq += zz * sin (qx * xx + qy * yy);

            }
	  //fq2_abs = fRq*fRq + fIq*fIq;
          fq2[count_q]+=fRq*fRq + fIq*fIq;
          count_q++;
        }
    }
  //return sum |fq|^2 
  return 0;
}





int calculate_fq2()
{
  double delta_q = 2.*M_PI/(double)NX;
  //double delta_q = 2.*M_PI/get_LX();
  
  double fRq=0, fIq=0, qx=0., qy=0., xx=0., yy=0., zz=0.;
  int count_q =0, nx=0, ny=0;
  //double sum_fq2=0, fq2_abs=0;
  double hcm=avg_hgt();
  double theta=0;
  //radius excluding R=0
  for (nx=1; nx<=NX; nx++)
    {
      
     
      theta=0;
      for (ny=0; ny<NY; ny++)
        {
	  theta = ny * 2.*M_PI/(double)NY;
          qx= (double) nx * delta_q * cos (theta);
          qy= (double) nx * delta_q * sin (theta);

          fRq=0.;
          fIq=0.;
          for(int i=0;i<N;i++)
            {
              xx = position[3*i];
              yy = position[3*i+1];
              zz = position[3*i+2]-hcm;


              fRq += zz * cos (qx * xx + qy * yy);
              fIq += zz * sin (qx * xx + qy * yy);

            }
	  
	  //average over N differrent directions
          fq2[count_q]+=(fRq*fRq + fIq*fIq)/(double)NY; //devide by N to avoid large number from sum 
        }
      count_q++;
    }
  return 0;
}




int avg_fq2(FILE *fq2data,int tot_frames)
{
  if (tot_frames == 0)
    {
      printf("Average fq2 computation is dividing by Zero\n");
      return 0;
    }
  double delta_q = 2.*M_PI/((double)NX);
  double qx=0., qy=0.; 
  int count_q =0, nx=0, ny=0, vol=(double)((NX)*(NY)*(NX)*(NY));
  for (nx=1; nx<=NX/2; nx++)
    {
      
      qx=delta_q * nx;
      
      
      fprintf(fq2data, "%.8f\t%.8f\n", qx, fq2[count_q]/(vol)/tot_frames);
      count_q++;
    }
  
  
   
   return 0;
}



int avg_fq2_1D(FILE *fq2data,int tot_frames)
{
   if (tot_frames == 0)
   {
        printf("Average fq2 computation is dividing by Zero\n");
        return 0;
   }
   double delta_q = 2.*M_PI/((double)NX);
   double qx=0., qy=0.;
   int count_q =0, nx=0, ny=0;
   for (nx=-NX/2; nx<NX/2; nx++)
    {
          qx=delta_q * (double) nx;

          fprintf(fq2data, "%.8f\t%.8f\n", qx, fq2[count_q]/(NX*NX)/tot_frames);
          count_q++;

    }

   return 0;
}



int calculate_fq2_symmetrized()
{
  //under development. save ~1/2 computing time 
  //double delta_q = 2.*M_PI/(double)NX;
  double delta_q = 2.*M_PI/get_LX();
  
  double fRq=0, fIq=0, qx=0., qy=0., xx=0., yy=0., zz=0.;
  int count_q =0, nx=0, ny=0, ii=0;
  
  //double sum_fq2=0, fq2_abs=0;
  //use symmetry to sample only 1/2 of whole q-space
  double hcm=avg_hgt();
  //first exclude the modes that include qx=0 or qy=0
  //use top half of the space 
  for (nx=-NX/2; nx<=NX/2; nx++)
    {
      for (ny=1; ny<=NY/2; ny++)
        {
	  if (nx!=0)
	    {
	      qx=delta_q * (double) nx;
	      qy=delta_q * (double) ny;
	      
	      fRq=0.;
	      fIq=0.;
	      for(int i=0;i<N;i++)
		{
		  xx = position[3*i];
		  yy = position[3*i+1];
		  zz = position[3*i+2]-hcm;
		  
		  
		  fRq += zz * cos (qx * xx + qy * yy);
		  fIq += zz * sin (qx * xx + qy * yy);
		  
		}
	      //fq2_abs = fRq*fRq + fIq*fIq;
	      fq2[count_q]+=(fRq*fRq + fIq*fIq);
	      count_q++;
	    }
	  //use symmetry 
	  
	      
	  
	  
	}
    }
  
  //next use symmetry to fill out f(-qx, -qy) using f(qx, qy) calculated ealier 
 
  for (int i=0; i<count_q; i++)
    {
      fq2[i+2*count_q]=fq2[i];
      
    }
  //update count_q as now we doubled our fq2
  count_q = 2*count_q; 

  //next calculate along qx and qy line 
  for (nx=-NX/2; nx<=NX/2; nx++)
    {
      ny=0; 
      if (nx!=0)
	{
	  qx=delta_q * (double) nx;
	  qy=delta_q * (double) ny;
      
	  fRq=0.;
	  fIq=0.;
	  for(int i=0;i<N;i++)
	    {
              xx = position[3*i];
              yy = position[3*i+1];
              zz = position[3*i+2]-hcm;
	      
	      
              fRq += zz * cos (qx * xx + qy * yy);
              fIq += zz * sin (qx * xx + qy * yy);
	      
	    }
	  //fq2_abs = fRq*fRq + fIq*fIq;
	  fq2[count_q]+=fRq*fRq + fIq*fIq;
	  count_q++;
	}
    }

   for (ny=-NY/2; ny<=NY/2; ny++)
    {
      nx=0; 
      if (ny!=0)
	{

	  qx=delta_q * (double) nx;
	  qy=delta_q * (double) ny;
	  
	  fRq=0.;
	  fIq=0.;
	  for(int i=0;i<N;i++)
	    {
              xx = position[3*i];
              yy = position[3*i+1];
              zz = position[3*i+2]-hcm;
	      
	      
              fRq += zz * cos (qx * xx + qy * yy);
              fIq += zz * sin (qx * xx + qy * yy);
	      
	    }
	  //fq2_abs = fRq*fRq + fIq*fIq;
	  fq2[count_q]+=fRq*fRq + fIq*fIq;
	  count_q++;
	}
    }
      
   //lastly qx=0, qy=0
   qx=0.;
   qy=0.; 
   fRq=0.;
   fIq=0.;
   for(int i=0;i<N;i++)
     {
       xx = position[3*i];
       yy = position[3*i+1];
       zz = position[3*i+2]-hcm;
       
       
       fRq += zz * cos (qx * xx + qy * yy);
       fIq += zz * sin (qx * xx + qy * yy);
       
     }
   //fq2_abs = fRq*fRq + fIq*fIq;
   fq2[count_q]+=fRq*fRq + fIq*fIq;
   count_q++; 
  return 0;
}





