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


double find_xmax()
{
  double xmax=-100000., xx=0.;
  for(int i=0;i<N;i++)
    {
      xx = position[3*i];
      if (xx>xmax)
        {
          xmax=xx;
        }

    }
  return (xmax);
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


double find_ymax()
{
  double ymax=-100000., yy=0;
  for(int i=0;i<N;i++)
    {
      yy = position[3*i+1];
      if(yy>ymax)
        {
          ymax=yy;

        }
    }
  return (ymax);
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


double build_fmesh()
{
  //(0,0) coordinates
  double LX_MIN=find_xmin(), LY_MIN=find_ymin(), LX_MAX=find_xmax(), LY_MAX=find_ymax(), hcm=avg_hgt();
  double delta =DELTA, xx=0, yy=0, zz=0; 
  int nx=0, ny=0, index=0; 
  int NXX=0, NYY=0;
  NXX = int((LX_MAX-LX_MIN)/delta);
  NYY = int((LY_MAX-LY_MIN)/delta);
  //printf("NXX, NYY: %d\t%d\n", NXX, NYY);
  /*
    grid NX * NY:
    
    
    
    NX NX+1 ...     2NX-1 
    0 1 2 3 4 5 ... NX-1
  */
  
   for(int i=0;i<(NXX*NYY); i++)
    { 
      //printf ("%d\n", i);
      fmesh[i]=0.;
      countmesh[i]=0;
      //                    
     }
      //
  
  
  
  for(int i=0; i<N; i++)
    {
      //relative position to (0,0)
      xx = position[3*i]-LX_MIN;
      yy = position[3*i+1]-LY_MIN;
      zz = position[3*i+2]-hcm;
      
      //convert to grid location
      nx = int(xx);
      ny = int(yy); 
      //if (nx>NX)
	//printf("%d\t%lf\n", nx, xx);
      //convert nx, ny to index in 1D array 
      index = nx%NXX+ny*NXX;
      //if (index>NX*NY)
      //printf ("%d\t%d\t%d\n", nx, ny, index); 
      fmesh[index]+=zz;
      countmesh[index]+=1;
    }
  
  //printf("%d", NX*NY);  
  //average atom height in each mesh 
  for(int i=0; i<NXX*NYY; i++)
    {
      if (countmesh[i]!=0)
	{
	    fmesh[i]=fmesh[i]/((double)countmesh[i]);
            //printf("%lf\n", fmesh[i]);
	}
       //printf("%lf\n", fmesh[i]); 
    }
  
  //printf("DONE"); 
  //return 0;


//double calculate_sum_dfdxi2()

  //calculate <(df/dxi)^2>
  double sum=0.;
  //delta_spacing grid = 1 unit 
  //excluding boundaries 
  
  for (nx=0; nx<NXX-1; nx++)
    {
      for (ny=0; ny<NYY-1; ny++)
        {
	  
	  index = nx + ny*NXX;
	  
	  //(df/dx)**2 + (df/dy)**2
	  sum  += (fmesh[index+1]-fmesh[index])*(fmesh[index+1]-fmesh[index]) + (fmesh[index+NXX]-fmesh[index])*(fmesh[index+NXX]-fmesh[index]);
	  
	}
    }

  //x-boundaries (last column)
  for (ny=0; ny<NYY-1; ny++)
    {
      //(df/dx)**2  + (df/dy)**2
      //last column (NX-1), use PBC for x 
      index = (NXX-1)+ny*NXX;
      sum  += (fmesh[index-NXX+1]-fmesh[index])*(fmesh[index-NXX+1]-fmesh[index])+(fmesh[index+NXX]-fmesh[index])*(fmesh[index+NXX]-fmesh[index]);
    }
  
  //y-boundaries
  for (nx=0; nx<NXX-1; nx++)
    {
      //(df/dx)**2 + (df/dy)**2
      //last row (NY-1), use PBC for y 
      index = nx+(NYY-1)*NXX;
      sum  += (fmesh[index+1]-fmesh[index])*(fmesh[index+1]-fmesh[index]) +(fmesh[index]-fmesh[nx])*(fmesh[index]-fmesh[nx]);
    }
  //corner
  sum  += (fmesh[NXX*NYY-NXX]-fmesh[NXX*NYY-1])*(fmesh[NXX*NYY-NXX]-fmesh[NXX*NYY-1]) +(fmesh[NXX-1]-fmesh[NXX*NYY-1])*(fmesh[NXX-1]-fmesh[NXX*NYY-1]);
  //count_n+=4;
  //printf("-------%d--------\n", count_n); 
  if (CALC==1)
    {
      return sum/(NXX)/(NYY)/delta/delta;
    }
  else
    { //CALC==2  -- use area T=0, A_0
      return sum/delta/delta;

    }
  
  //return sum/NXX/NYY;
}

double build_fmesh_cd4()
{
  //(0,0) coordinates
  double LX_MIN=find_xmin(), LY_MIN=find_ymin(), LX_MAX=find_xmax(), LY_MAX=find_ymax(), hcm=avg_hgt();
  double delta = DELTA, xx=0, yy=0, zz=0, diff=0; 
  int nx=0, ny=0, index=0, nxx=0, nyy=0; 
  int NXX=0, NYY=0;
  //delta = (LX_MAX-LX_MIN)/((double)NX);
  NXX = int(ceil((LX_MAX-LX_MIN)/delta))+1;
  NYY = int(ceil((LY_MAX-LY_MIN)/delta))+1;
  
  //printf("NXX, NYY: %d\t%d\n", NXX, NYY);
  /*
    grid NX * NY:
    
    
    
    NX NX+1 ...     2NX-1 
    0 1 2 3 4 5 ... NX-1
  */
  
  for(int i=0;i<(NXX*NYY); i++)
    { 
      //printf ("%d\n", i);
      fmesh[i]=0.;
      countmesh[i]=0;
      //                    
    }
  //
  
  
  
  for(int i=0; i<N; i++)
    {
      //relative position to (0,0)
      xx = position[3*i]-LX_MIN + 0*delta/2.;
      yy = position[3*i+1]-LY_MIN + 0*delta/2.;
      zz = position[3*i+2]-hcm;
      
      //convert to grid location
      nx = int(xx/delta);
      ny = int(yy/delta); 
      //if (nx>NX)
      //printf("%d\t%lf\n", nx, xx);
      //convert nx, ny to index in 1D array 
      index = nx%NXX+ny*NXX;
      //if (index>NX*NY)
      //printf ("%d\t%d\t%d\n", nx, ny, index); 
      
      if (EXCL==0) //include all, exclude dilations==NO 
	{
	  fmesh[index]+=zz;
	  countmesh[index]+=1;
	}
      else if (EXCL==1) //exlude dilations 
	{
	  if (i%2==0 and (i/NY)%2==0)
	    {
	      //print buckled site for debugging 
	      //printf("id: %d\t", i);
	      //use mid point (forward for inner and backward for boundaries)
	      
	      if (nx<NXX-1 and ny<NYY-1)
		{
		  fmesh[index]+=(position[3*(i+1)+2]+position[3*(i+NY)+2])/2-hcm;
		  countmesh[index]+=1;
		}
	      else
		{
		  fmesh[index]+=(position[3*(i-1)+2]+position[3*(i-NY)+2])/2-hcm;
		  countmesh[index]+=1;
		}
	      
	    }
	  else
	    {
	      fmesh[index]+=zz;
	      countmesh[index]+=1;
	    }
	}

    }
  
  //printf("%d", NX*NY);  
  //average atom height in each mesh 
  for(int i=0; i<NXX*NYY; i++)
    {
      if (countmesh[i]!=0)
	{
	  fmesh[i]=fmesh[i]/((double)countmesh[i]);
	  //printf("%lf\n", fmesh[i]);
	}
      else
	{
	  fmesh[i]=0.;
	  
	}
      //printf("%lf\n", fmesh[i]); 
    }
  
  //printf("DONE"); 
  //return 0;
  
  
  //double calculate_sum_dfdxi2()

  //calculate <(df/dxi)^2>
  double sum=0., sum2=0;
  int count_n=0; 
  //delta_spacing grid = 1 unit 
  //excluding boundaries 
  
  for (nx=2; nx<NXX-2; nx++)
    {
      for (ny=2; ny<NYY-2; ny++)
        {
	  
	  index = nx + ny*NXX;
	  
	  //(df/dx)**2 + (df/dy)**2
	  diff  =-fmesh[index+2]+8*fmesh[index+1]-8*fmesh[index-1]+fmesh[index-2];
	  // (fmesh[index+1]-fmesh[index-1])*(fmesh[index+1]-fmesh[index-1]) + (fmesh[index+NXX]-fmesh[index-NXX])*(fmesh[index+NXX]-fmesh[index-NXX]);
	  sum += diff*diff/12./12./delta/delta;
	  diff = -fmesh[index+2*NXX]+8*fmesh[index+NXX]-8*fmesh[index-NXX]+fmesh[index-2*NXX];;
	  sum += diff*diff/12./12./delta/delta;
        //count_n++;	  
	}
    }
  //inner boundaries 
  //nx=1;
  //ny=1;

  //index = nx + ny*NXX;
  
  //sum+= ((fmesh[index+1]-fmesh[index-1])*(fmesh[index+1]-fmesh[index-1]) + (fmesh[index+NXX]-fmesh[index-NXX])*(fmesh[index+NXX]-fmesh[index-NXX]))/4./delta/delta;


  
  
  //nx=NXX-1;
  //ny=NYY-1;
  
  //sum+= ((fmesh[index+1]-fmesh[index-1])*(fmesh[index+1]-fmesh[index-1]) + (fmesh[index+NXX]-fmesh[index-NXX])*(fmesh[index+NXX]-fmesh[index-NXX]))/4./delta/delta;
  
  



  
  
  
  if (PBC==1)
    {
      //x-boundaries (right column)
      for (ny=1; ny<NYY-1; ny++)
	{
	  //(df/dx)**2  + (df/dy)**2
	  //right column (NX-1), use PBC for x 
	  index = (NXX-1)+ny*NXX;
	  sum2  += (fmesh[index-NXX+1]-fmesh[index-1])*(fmesh[index-NXX+1]-fmesh[index-1])+(fmesh[index+NXX]-fmesh[index-NXX])*(fmesh[index+NXX]-fmesh[index-NXX]);
	  
	  //left column
	  index = ny*NXX;
	  sum2  += (fmesh[index+1]-fmesh[index+NXX-1])*(fmesh[index+1]-fmesh[index+NXX-1])+(fmesh[index+NXX]-fmesh[index-NXX])*(fmesh[index+NXX]-fmesh[index-NXX]);
	  
	  //count_n +=2;
	}
      
      
      //y-boundaries
      for (nx=1; nx<NXX-1; nx++)
	{
	  //(df/dx)**2 + (df/dy)**2
	  //top row (NY-1), use PBC for y 
	  index = nx+(NYY-1)*NXX;
	  sum2  += (fmesh[index+1]-fmesh[index-1])*(fmesh[index+1]-fmesh[index-1]) +(fmesh[nx]-fmesh[index-NXX])*(fmesh[nx]-fmesh[index-NXX]);
	  
	  //bottom row
	  index = nx;
	  sum2  += (fmesh[index+1]-fmesh[index-1])*(fmesh[index+1]-fmesh[index-1]) +(fmesh[nx+NXX]-fmesh[nx+(NYY-1)*NXX])*(fmesh[nx+NXX]-fmesh[nx+(NYY-1)*NXX]);
	  //count_n+=2; 
	}
      //corner
      //top right
      sum2  += (fmesh[NXX*NYY-NXX]-fmesh[NXX*NYY-2])*(fmesh[NXX*NYY-NXX]-fmesh[NXX*NYY-2]) +(fmesh[NXX-1]-fmesh[NXX*NYY-1-NXX])*(fmesh[NXX-1]-fmesh[NXX*NYY-1-NXX]);
      //top left
      sum2  += (fmesh[NXX*NYY-NXX+1]-fmesh[NXX*NYY-1])*(fmesh[NXX*NYY-NXX+1]-fmesh[NXX*NYY-1]) +(fmesh[0]-fmesh[NXX*(NYY-1)-NXX])*(fmesh[0]-fmesh[NXX*(NYY-1)-NXX]);
      //bottom left 
      sum2  += (fmesh[1]-fmesh[NXX-1])*(fmesh[1]-fmesh[NXX-1]) +(fmesh[NXX]-fmesh[NXX*(NYY-1)])*(fmesh[NXX]-fmesh[NXX*(NYY-1)]);
      
      sum2  += (fmesh[0]-fmesh[NXX-2])*(fmesh[0]-fmesh[NXX-2]) +(fmesh[NXX-1+NXX]-fmesh[NXX*NYY-1])*(fmesh[NXX-1+NXX]-fmesh[NXX*NYY-1]);
      
    }
  sum2=sum2/4./delta/delta;
  //count_n+=4;
  //printf("-------%d--------\n", count_n); 
  if (CALC==1)
    {
      return (sum+sum2)/(NXX)/(NYY);
    }
  else
    { //CALC==2  -- use area T=0, A_0
      return sum+sum2;
      
    }
}


double build_fmesh_cd()
{
  //(0,0) coordinates
  double LX_MIN=find_xmin(), LY_MIN=find_ymin(), LX_MAX=find_xmax(), LY_MAX=find_ymax(), hcm=avg_hgt();
  double delta = DELTA, xx=0, yy=0, zz=0; 
  int nx=0, ny=0, index=0, nxx=0, nyy=0; 
  int NXX=0, NYY=0;
  
  if (delta<0){
  	delta = (LX_MAX-LX_MIN)/((double)NX);
  	NXX = NX; //int(ceil((LX_MAX-LX_MIN)/delta))+1;
  	NYY = NY; //int(ceil((LY_MAX-LY_MIN)/delta))+1;
  }
  else{
  NXX = int(ceil((LX_MAX-LX_MIN)/delta))+1;
  NYY = int(ceil((LY_MAX-LY_MIN)/delta))+1;
  }
    
  //printf("NXX, NYY: %d\t%d\n", NXX, NYY);
  /*
    grid NX * NY:
    
    
    
    NX NX+1 ...     2NX-1 
    0 1 2 3 4 5 ... NX-1
  */
  
   for(int i=0;i<(NXX*NYY); i++)
    { 
      //printf ("%d\n", i);
      fmesh[i]=0.;
      countmesh[i]=0;
      //                    
     }
      //
  
  
  
  for(int i=0; i<N; i++)
    {
      //relative position to (0,0)
      xx = position[3*i]-LX_MIN + 0*delta/2.;
      yy = position[3*i+1]-LY_MIN + 0*delta/2.;
      zz = position[3*i+2]-hcm;
      
      //convert to grid location
      nx = int(xx/delta);
      ny = int(yy/delta); 
      //if (nx>NX)
	//printf("%d\t%lf\n", nx, xx);
      //convert nx, ny to index in 1D array 
      index = nx%NXX+ny*NXX;
      //if (index>NX*NY)
      //printf ("%d\t%d\t%d\n", nx, ny, index); 
      
      if (EXCL==0) //include all, exclude dilations==NO 
	{
	  fmesh[index]+=zz;
	  countmesh[index]+=1;
	}
      else if (EXCL==1) //exlude dilations 
	{
	  if (i%2==0 and (i/NY)%2==0)
	    {
	      //print buckled site for debugging 
	      //printf("id: %d\t", i);
	      //use mid point (forward for inner and backward for boundaries)
	      
	      if (nx<NXX-1 and ny<NYY-1)
		{
		  fmesh[index]+=(position[3*(i+1)+2]+position[3*(i+NY)+2])/2-hcm;
		  countmesh[index]+=1;
		}
	      else
		{
		  fmesh[index]+=(position[3*(i-1)+2]+position[3*(i-NY)+2])/2-hcm;
		  countmesh[index]+=1;
		}
	      
	    }
	  else
	    {
	      fmesh[index]+=zz;
	      countmesh[index]+=1;
	    }
	}

    }
  
  //printf("%d", NX*NY);  
  //average atom height in each mesh 
  for(int i=0; i<NXX*NYY; i++)
    {
      if (countmesh[i]!=0)
	{
	  fmesh[i]=fmesh[i]/((double)countmesh[i]);
	  //printf("%d\t%lf\t%d\n", i, fmesh[i], NXX);
	}
      else
	{
	  fmesh[i]=0.; //(fmesh[i-1]+fmesh[i+1])/2.;
	  
	}
       //printf("%lf\n", fmesh[i]); 
    }
  
  //printf("DONE"); 
  //return 0;


//double calculate_sum_dfdxi2()

  //calculate <(df/dxi)^2>
  double sum=0.;
  int count_n=0; 
  //delta_spacing grid = 1 unit 
  //excluding boundaries 
  
  for (nx=1; nx<NXX-1; nx++)
    {
      for (ny=1; ny<NYY-1; ny++)
        {
	  
	  index = nx + ny*NXX;
	  
	  //(df/dx)**2 + (df/dy)**2
	  sum  += (fmesh[index+1]-fmesh[index-1])*(fmesh[index+1]-fmesh[index-1]) + (fmesh[index+NXX]-fmesh[index-NXX])*(fmesh[index+NXX]-fmesh[index-NXX]);
        //count_n++;	  
	}
    }
  if (PBC==1)
    {
      //x-boundaries (right column)
      for (ny=1; ny<NYY-1; ny++)
	{
	  //(df/dx)**2  + (df/dy)**2
	  //right column (NX-1), use PBC for x 
	  index = (NXX-1)+ny*NXX;
	  sum  += (fmesh[index-NXX+1]-fmesh[index-1])*(fmesh[index-NXX+1]-fmesh[index-1])+(fmesh[index+NXX]-fmesh[index-NXX])*(fmesh[index+NXX]-fmesh[index-NXX]);
	  
	  //left column
	  index = ny*NXX;
	  sum  += (fmesh[index+1]-fmesh[index+NXX-1])*(fmesh[index+1]-fmesh[index+NXX-1])+(fmesh[index+NXX]-fmesh[index-NXX])*(fmesh[index+NXX]-fmesh[index-NXX]);
	  
	  //count_n +=2;
	}
      
      
  //y-boundaries
      for (nx=1; nx<NXX-1; nx++)
	{
	  //(df/dx)**2 + (df/dy)**2
	  //top row (NY-1), use PBC for y 
	  index = nx+(NYY-1)*NXX;
	  sum  += (fmesh[index+1]-fmesh[index-1])*(fmesh[index+1]-fmesh[index-1]) +(fmesh[nx]-fmesh[index-NXX])*(fmesh[nx]-fmesh[index-NXX]);
	  
	  //bottom row
	  index = nx;
	  sum  += (fmesh[index+1]-fmesh[index-1])*(fmesh[index+1]-fmesh[index-1]) +(fmesh[nx+NXX]-fmesh[nx+(NYY-1)*NXX])*(fmesh[nx+NXX]-fmesh[nx+(NYY-1)*NXX]);
	  //count_n+=2; 
	}
      //corner
      //top right
      sum  += (fmesh[NXX*NYY-NXX]-fmesh[NXX*NYY-2])*(fmesh[NXX*NYY-NXX]-fmesh[NXX*NYY-2]) +(fmesh[NXX-1]-fmesh[NXX*NYY-1-NXX])*(fmesh[NXX-1]-fmesh[NXX*NYY-1-NXX]);
      //top left
      sum  += (fmesh[NXX*NYY-NXX+1]-fmesh[NXX*NYY-1])*(fmesh[NXX*NYY-NXX+1]-fmesh[NXX*NYY-1]) +(fmesh[0]-fmesh[NXX*(NYY-1)-NXX])*(fmesh[0]-fmesh[NXX*(NYY-1)-NXX]);
      //bottom left 
      sum  += (fmesh[1]-fmesh[NXX-1])*(fmesh[1]-fmesh[NXX-1]) +(fmesh[NXX]-fmesh[NXX*(NYY-1)])*(fmesh[NXX]-fmesh[NXX*(NYY-1)]);
      
      sum  += (fmesh[0]-fmesh[NXX-2])*(fmesh[0]-fmesh[NXX-2]) +(fmesh[NXX-1+NXX]-fmesh[NXX*NYY-1])*(fmesh[NXX-1+NXX]-fmesh[NXX*NYY-1]);
  
    }
  //count_n+=4;
  //printf("-------%d--------\n", count_n); 
  if (CALC==1)
    {
      return sum/(NXX)/(NYY)/4/delta/delta;
    }
  else
    { //CALC==2  -- use area T=0, A_0
      return sum/4/delta/delta;

    }
}





double build_fmesh_cd_nodes()
{
  //(0,0) coordinates
  double LX_MIN=find_xmin(), LY_MIN=find_ymin(), LX_MAX=find_xmax(), LY_MAX=find_ymax(), hcm=avg_hgt();
  double delta = DELTA, xx=0, yy=0, zz=0; 
  int nx=0, ny=0, index=0; 
  int NXX=0, NYY=0;
  //delta = (LX_MAX-LX_MIN)/(double)NX;
  NXX = NX; //int((LX_MAX-LX_MIN)/delta);
  NYY = NY; //int((LY_MAX-LY_MIN)/delta);
  
  //printf("NXX, NYY: %d\t%d\n", NXX, NYY);
  /*
    grid NX * NY:
    
    
    
    NX NX+1 ...     2NX-1 
    0 1 2 3 4 5 ... NX-1
  */
  
  //		for(int i=0;i<(NXX*NYY); i++)
   
  double sum=0.;
  int count_n=0; 
  //delta_spacing grid = 1 unit 
  //excluding boundaries 
  
  for (nx=1; nx<NXX-1; nx++)
    {
      for (ny=1; ny<NYY-1; ny++)
        {
	  
	  index = nx + ny*NXX;
	  
	  
	  //xx = position[3*i];
	  //(df/dx)**2 + (df/dy)**2
	  sum  += (position[3*(index+1)+2]-position[3*(index-1)+2])*(position[3*(index+1)+2]-position[3*(index-1)+2])/(position[3*(index+1)]-position[3*(index-1)])/(position[3*(index+1)]-position[3*(index-1)]);

///(position[3*(index)]-position[3*(index-1)]);
	  
	  sum  += (position[3*(index+NXX)+2]-position[3*(index-NXX)+2])*(position[3*(index+NXX)+2]-position[3*(index-NXX)+2])/(position[3*(index+NXX)+1]-position[3*(index-NXX)+1])*(position[3*(index+NXX)+1]-position[3*(index-NXX)+1]);

//(position[3*index+1]-position[3*(index-NXX)+1]);
	  //count_n++;	  
	  
	}
    }
  
  //x-boundaries (right column)
  //for (ny=1; ny<NYY-1; ny++)
  //  {
  //(df/dx)**2  + (df/dy)**2
  //right column (NX-1), use PBC for x 
  //  index = (NXX-1)+ny*NXX;
  //  sum  += (fmesh[index-NXX+1]-fmesh[index-1])*(fmesh[index-NXX+1]-fmesh[index-1])+(fmesh[index+NXX]-fmesh[index-NXX])*(fmesh[index+NXX]-fmesh[index-NXX]);

  //left column
  //  index = ny*NXX;
  //  sum  += (fmesh[index+1]-fmesh[index+NXX-1])*(fmesh[index+1]-fmesh[index+NXX-1])+(fmesh[index+NXX]-fmesh[index-NXX])*(fmesh[index+NXX]-fmesh[index-NXX]);

  //count_n +=2;
  //}
  
  
  //y-boundaries
  //for (nx=1; nx<NXX-1; nx++)
  // {
  //(df/dx)**2 + (df/dy)**2
  //top row (NY-1), use PBC for y 
  //  index = nx+(NYY-1)*NXX;
  //  sum  += (fmesh[index+1]-fmesh[index-1])*(fmesh[index+1]-fmesh[index-1]) +(fmesh[nx]-fmesh[index-NXX])*(fmesh[nx]-fmesh[index-NXX]);
  
  //bottom row
  //  index = nx;
  //sum  += (fmesh[index+1]-fmesh[index-1])*(fmesh[index+1]-fmesh[index-1]) +(fmesh[nx+NXX]-fmesh[nx+(NYY-1)*NXX])*(fmesh[nx+NXX]-fme//sh[nx+(NYY-1)*NXX]);
  //count_n+=2; 
  //}
  //corner
  //top right
  //sum  += (fmesh[NXX*NYY-NXX]-fmesh[NXX*NYY-2])*(fmesh[NXX*NYY-NXX]-fmesh[NXX*NYY-2]) +(fmesh[NXX-1]-fmesh[NXX*NYY-1-NXX])*(fmesh[NXX-1]-fmesh[NXX*NYY-1-NXX]);
  //top left
  //sum  += (fmesh[NXX*NYY-NXX+1]-fmesh[NXX*NYY-1])*(fmesh[NXX*NYY-NXX+1]-fmesh[NXX*NYY-1]) +(fmesh[0]-fmesh[NXX*(NYY-1)-NXX])*(fmesh[0]-fmesh[NXX*(NYY-1)-NXX]);
  //bottom left 
  //sum  += (fmesh[1]-fmesh[NXX-1])*(fmesh[1]-fmesh[NXX-1]) +(fmesh[NXX]-fmesh[NXX*(NYY-1)])*(fmesh[NXX]-fmesh[NXX*(NYY-1)]);
  
  //sum  += (fmesh[0]-fmesh[NXX-2])*(fmesh[0]-fmesh[NXX-2]) +(fmesh[NXX-1+NXX]-fmesh[NXX*NYY-1])*(fmesh[NXX-1+NXX]-fmesh[NXX*NYY-1]);
  //count_n+=4;
  //printf("-------%d--------\n", count_n); 
  if (CALC==3)
    {
      return sum/(NXX)/(NYY);
    }
  else if (CALC==4)
    { //CALC==2  -- use area T=0, A_0
      return sum;
      
    }
}




//to calculate <f^2> for debugging purposes 
double build_fmesh_f2()
{
  //(0,0) coordinates
  double LX_MIN=find_xmin(), LY_MIN=find_ymin(), LX_MAX=find_xmax(), LY_MAX=find_ymax(), hcm=avg_hgt();
  double delta = 1., xx=0, yy=0, zz=0; 
  int nx=0, ny=0, index=0; 
  int NXX=0, NYY=0;
  NXX = int((LX_MAX-LX_MIN)/delta);
  NYY = int((LY_MAX-LY_MIN)/delta);
  //printf("%d\t%d\n", NXX, NYY);
  /*
    grid NX * NY:
    
    
    
    NX NX+1 ...     2NX-1 
    0 1 2 3 4 5 ... NX-1
  */
  
   for(int i=0;i<(NXX*NYY); i++)
    { 
      //printf ("%d\n", i);
      fmesh[i]=0.;
      countmesh[i]=0;
      //                    
     }
      //
  
  
  
  for(int i=0; i<N; i++)
    {
      //relative position to (0,0)
      xx = position[3*i]-LX_MIN;
      yy = position[3*i+1]-LY_MIN;
      zz = position[3*i+2]-hcm;
      
      //convert to grid location
      nx = int(xx);
      ny = int(yy); 
      //if (nx>NX)
	//printf("%d\t%lf\n", nx, xx);
      //convert nx, ny to index in 1D array 
      index = nx%NXX+ny*NXX;
      //if (index>NX*NY)
      //printf ("%d\t%d\t%d\n", nx, ny, index); 
      fmesh[index]+=zz;
      countmesh[index]+=1;
    }
  
  //printf("%d", NX*NY);  
  //average atom height in each mesh 
  for(int i=0; i<NXX*NYY; i++)
    {
      if (countmesh[i]!=0)
	{
	    fmesh[i]=fmesh[i]/((double)countmesh[i]);
            //printf("%lf\n", fmesh[i]);
	}
       //printf("%lf\n", fmesh[i]); 
    }
  
  //printf("DONE"); 
  //return 0;


//double calculate_sum_dfdxi2()

  //calculate <(df/dxi)^2>
  double sum=0.;
  //int nx=0, ny=0, index=0;
  //delta_spacing grid = 1 unit 
  //excluding boundaries 
  
  for (nx=0; nx<NXX-1; nx++)
    {
      for (ny=0; ny<NYY-1; ny++)
        {
	  
	  index = nx + ny*NXX;
	  
	  //(df/dx)**2 + (df/dy)**2
	  sum +=fmesh[index]*fmesh[index];
	  //sum  += (fmesh[index+1]-fmesh[index])*(fmesh[index+1]-fmesh[index]) + (fmesh[index+NXX]-fmesh[index])*(fmesh[index+NXX]-fmesh[index]);
	  
	}
    }

  //x-boundaries (last column)
  for (ny=0; ny<NYY-1; ny++)
    {
      //(df/dx)**2  + (df/dy)**2
      //last column (NX-1), use PBC for x 
      index = (NXX-1)+ny*NXX;
      sum += fmesh[index]*fmesh[index];
      //sum  += (fmesh[index-NXX+1]-fmesh[index])*(fmesh[index-NXX+1]-fmesh[index])+(fmesh[index+NXX]-fmesh[index])*(fmesh[index+NXX]-fmesh[index]);
    }
  
  //y-boundaries
  for (nx=0; nx<NXX-1; nx++)
    {
      //(df/dx)**2 + (df/dy)**2
      //last row (NY-1), use PBC for y 
      index = nx+(NYY-1)*NXX;
      sum += fmesh[index]*fmesh[index];
      //sum  += (fmesh[index+1]-fmesh[index])*(fmesh[index+1]-fmesh[index]) +(fmesh[index]-fmesh[nx])*(fmesh[index]-fmesh[nx]);
    }
  //corner
  sum +=fmesh[NXX*NYY-1]*fmesh[NXX*NYY-1];
  //sum  += (fmesh[NXX*NYY-NXX]-fmesh[NXX*NYY-1])*(fmesh[NXX*NYY-NXX]-fmesh[NXX*NYY-1]) +(fmesh[NXX-1]-fmesh[NXX*NYY-1])*(fmesh[NXX-1]-fmesh[NXX*NYY-1]);
  return sum/NXX/NYY;
}

