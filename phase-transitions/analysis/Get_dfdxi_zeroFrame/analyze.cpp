#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "analyze.h"

double fmesh[NMAX];
int countmesh[NMAX];
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


double find_xmin()
{
  double xmin=100000., xx=0.;
  for(int i=0;i<N;i++)
    {
      xx = position[3*i];
      if (xx<xmin)
	{
	  xmin=xx;
	}
      
    }
  return (xmin);
}


double find_ymin()
{
  double ymin=100000., yy=0;
  for(int i=0;i<N;i++)
    {
      yy = position[3*i+1];
      if(yy<ymin)
	{
	  ymin=yy;
	  
	}
    }
  return (ymin);
}



int initialize_fmesh()
{
  for(int i=0;i<(NX*NY); i++)
    { 
      //printf ("%d\n", i);
      fmesh[i]=0.;
      countmesh[i]=0;
	
    }
  return 0;
}


int build_fmesh()
{
  //(0,0) coordinates
  double LX=find_xmin(), LY=find_ymin(), hcm=avg_hgt();
  double delta = 1., xx=0, yy=0, zz=0; 
  int nx=0, ny=0, index=0; 
  /*
    grid NX * NY:
    
    
    
    NX NX+1 ...     2NX-1 
    0 1 2 3 4 5 ... NX-1
  */
  
  
  
  
  
  for(int i=0; i<N; i++)
    {
      //relative position to (0,0)
      xx = position[3*i]-LX;
      yy = position[3*i+1]-LY;
      zz = position[3*i+2]-hcm;
      
      //convert to grid location
      nx = int(xx);
      ny = int(yy); 
      if (nx>NX)
	printf("%d\t%lf\n", nx, xx);
      //convert nx, ny to index in 1D array 
      index = nx%NX+ny*NX;
      if (index>NX*NY)
      printf ("%d\t%d\t%d\n", nx, ny, index); 
      fmesh[index]+=zz;
      countmesh[index]+=1;
    }
  
  //printf("%d", NX*NY);  
  //average atom height in each mesh 
  for(int i=0; i<NX*NY; i++)
    {
      if (countmesh[i]!=0)
	{
	    fmesh[i]=fmesh[i]/((double)countmesh[i]);
            //printf("%lf\n", fmesh[i]);
	}
       //printf("%lf\n", i, fmesh[i]); 
    }
  
  //printf("DONE"); 
  return 0;
}


double calculate_sum_dfdxi2()
{
  //calculate <(df/dxi)^2>
  double sum=0.;
  int nx=0, ny=0, index=0;
  //delta_spacing grid = 1 unit 
  //excluding boundaries 
  
  for (nx=0; nx<NX-1; nx++)
    {
      for (ny=0; ny<NY-1; ny++)
        {
	  
	  index = nx + ny*NX;
	  
	  //(df/dx)**2 + (df/dy)**2
	  
	  sum  += (fmesh[index+1]-fmesh[index])*(fmesh[index+1]-fmesh[index]) + (fmesh[index+NX]-fmesh[index])*(fmesh[index+NX]-fmesh[index]);
	  
	}
    }

  //x-boundaries
  for (ny=0; ny<NY; ny++)
    {
      //(df/dx)**2 ONLY in x-direction only
      //last column (NX-1)
      index = (NX-1)+ny*NX;
      sum  += (fmesh[index-NX]-fmesh[index])*(fmesh[index-NX]-fmesh[index]);
    }
  
  //y-boundaries
  for (nx=0; nx<NX; nx++)
    {
      //(df/dx)**2 ONLY in y-direction only
      //last row (NY-1)
      index = nx+(NY-1)*NX;
      sum  += (fmesh[index]-fmesh[nx])*(fmesh[index]-fmesh[nx]);
    }
  
  return sum;
}



