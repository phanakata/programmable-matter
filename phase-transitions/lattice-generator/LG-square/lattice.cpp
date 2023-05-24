#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <limits>
#include <string.h>
#include <iostream>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>
#include "lattice_variables.h"
#include "lattice.h"
#include "stdint.h"


int bond_mat[NMAX][NMAX];
int num_bonds = 0;
int dihedrals[MAXDIHEDRALS][4];
int cnt_dihedrals=0;
int particle_id[NMAX];
int bondType;
/*--------------------------------------------------------------------------------------*/
/*      Function for setting up Initial lattice coordinates                             */
/*--------------------------------------------------------------------------------------*/
int initialLatticeStruct(latticeStruct *ip,  int size)
{
  int k = 0;/*keeping track of coords along x-axis*/
  int j; 
  for (int i = 0; i < size; i++)
  {
        //pzh need to shift x and y to fit all atoms in the box
	j = i/NX;
        ip->x[i] = (double)(k * a-(double)NX/2);
        k++;
        k = (k%NX == 0) ? 0 : k; /*resetting k to 0 once the lattice has been traversed in the x direction*/
        ip->y[i] = a  * (j-(double)NY/2);
        ip->z[i]=0; /*  Starting with flat configuration        */
        //ip->z[NX/2]=grandom(0,ZRAND);//Moving one node randomly
   }

  if (defect==2)
    {
      
      for(int i=edgex; i<NX-edgex;i+=spx)
        {       
	  for (int j=edgey; j<NY-edgey; j+=spy)
	    {
	      ip->z[j*NX+i]=0.1*pow(-1, (i-edgex)/spx+(j-edgey)/spy);
	    }     //bumps
	}
    }

// for negative defects give bumps to atoms between 4 defects 
 if (defect==-2)
    {

      for(int i=edgex+spx/2; i<NX-edgex;i+=spx)
        {
          for (int j=edgey+spy/2; j<NY-edgey; j+=spy)
            {
              ip->z[j*NX+i]=0.1*pow(-1, (i-(edgex+spx/2))/spx+(j-(edgey+spy/2))/spy);
            }     //bumps
        }
    }




  if (defect==4)
	{
	for (int i=0; i<LEN; i++)
		if (i==NX*(NY/2)+NX/2)
			{ip->z[i]=0.1;} 
	}
       // for(int i=3; i<NX-3;i+=spx)
        //{       for (int j=3; j<NY-3; j+=spy)
         //                       {ip->z[j*NX+i]=0.1*pow(-1, (i-3)/4+(j-3)/4);}     //bumps
       // }


	
   return 0;
}

/*----------------------------------------------------------------------------------------*/
/*Function for generating the Bond Matrix,Matrix[i][j] = 1,Bond between particles i and j */
/*----------------------------------------------------------------------------------------*/
int lattice_connectivity()
{

  /* Initializing elements of the Bond Matrix to 0	*/
  for(int i=0;i<LEN;i++)
  {
	for(int j=0;j<LEN;j++)
	{
		bond_mat[i][j]=0;
	}
  }
		
  for(int i=0;i<LEN;i++)
  {	
	/* LEFT or RIGHT Boundary */
	if (i%NX==0 or i%NX==NX-1)
	{
		//LEFT BOUNDARY 
		if (i%NX==0)
		{	
			//EVEN ROW 
			if ((i/NX)%2==0) 
			{
				if (i+NX<LEN)
                                	bond_mat[i][i+NX]=1;
				if (i+NX+1<LEN)
                                	bond_mat[i][i+NX+1]=2;
                        	if(i-NX >= 0)
                                	bond_mat[i][i-NX]=1;
				if(i-NX+1 >= 0)
                                	bond_mat[i][i-NX+1]=2;
                        	bond_mat[i][i+1]=1;
			}
			//ODD ROW
			else
			{
				if (i+NX<LEN)
					bond_mat[i][i+NX]=1;
				if(i-NX >= 0)
              				bond_mat[i][i-NX]=1;
           			bond_mat[i][i+1]=1;
			}
		}
		//RIGHT BOUNDARY 
		else
		{	
			//EVEN ROW
			if ((i/NX)%2==0)
                	{
                        	if (i+NX<LEN)
                                	bond_mat[i][i+NX]=1;
                        	if (i+NX-1<LEN)
                                	bond_mat[i][i+NX-1]=2;
                        	if(i-NX >= 0)
                                	bond_mat[i][i-NX]=1;
                        	if(i-NX-1 >= 0)
                                	bond_mat[i][i-NX-1]=2;
                        	bond_mat[i][i-1]=1;
                	}
			//ODD ROW
			else
			{	if (i+NX<LEN)
                                	bond_mat[i][i+NX]=1;
                        	if(i-NX >= 0)
                                	bond_mat[i][i-NX]=1;
                        	bond_mat[i][i-1]=1;
			}	
			
		}
	}
	/* INTERMEDIATE */
	else
	{
		/* ODD ROW and ODD COLUMN or EVEN ROW and EVEN COLUMN */
		if (((i/NX)%2==0 and (i%NX)%2==0) or ((i/NX)%2==1 and (i%NX)%2==1))
		{	
			bond_mat[i][i-1]=1;
                        bond_mat[i][i+1]=1;
			if(i+NX<LEN)
                                bond_mat[i][i+NX]=1;
			if(i+NX+1<LEN)
                                bond_mat[i][i+NX+1]=2;
			if(i+NX-1<LEN)
                                bond_mat[i][i+NX-1]=2;
                        if(i-NX>=0)
                                bond_mat[i][i-NX]=1;
			if(i-NX+1>=0)
                                bond_mat[i][i-NX+1]=2;
			if(i-NX-1>=0)
                                bond_mat[i][i-NX-1]=2;


			

		}
		/* EVEN ROW and ODD COLUMN or ODD ROW and EVEN COLUMN */
		else
		{
                	
			bond_mat[i][i-1]=1;
                        bond_mat[i][i+1]=1;
			if(i+NX<LEN)
                                bond_mat[i][i+NX]=1;
                        if(i-NX>=0)
                                bond_mat[i][i-NX]=1;

                }
	
	}
	
   }
   return 0;
}


int lattice_connectivity_pbc()
{
  int image_top, image_bottom;
  /* Initializing elements of the Bond Matrix to 0	*/
  for(int i=0;i<LEN;i++)
  {
	for(int j=0;j<LEN;j++)
	{
		bond_mat[i][j]=0;
	}
  }
		
  for(int i=0;i<LEN;i++)
    {
      /* 
	 x: left/right boundary 
	 y: top/bottom boundary
	 *: intermediate 
	 
	 x * * * * * x
	 x * * * * * x
	 x * * * * * x 
	 x * * * * * x
	 x * * * * * x
	 

       */
      
      /* LEFT or RIGHT Boundary */
      //need to debug this conditional for right boundary 
      if (i%NX==0 or ((i+1)%NX==0))
	{
	  //LEFT BOUNDARY 
	  if (i%NX==0)
	    {	
	      //EVEN ROW 
	      if ((i/NX)%2==0) 
		{
		  if (i!=0 and i!=NX*(NY-1))
		    {
		      //excluding corner so we gurantee i+NX<LEN and i-NX>=0
		      //if (i+NX<LEN)
		      
		      bond_mat[i][i+NX]=1;
		      bond_mat[i][i+NX+1]=2;
		      bond_mat[i][i+1]=1;
		      bond_mat[i][i-NX+1]=2;
		      bond_mat[i][i-NX]=1;
		      //PBC image
		      bond_mat[i][i-1]=2;
		      bond_mat[i][i-1+NX]=1;
		      bond_mat[i][i-1+NX+NX]=2;
		      

		    }
		  else
		    {
		      //no top right (left?) corner for even ROW 
		      //bottom left corner 
		      bond_mat[i][i+NX]=1; 
		      bond_mat[i][i+NX+1]=2;
		      bond_mat[i][i+1]=1;
		      bond_mat[i][i+NX*(NY-1)+1]=2; 
		      bond_mat[i][i+NX*(NY-1)]=1;
		      bond_mat[i][i+NX*NY-1]=2;// 0 and 47 - 2 
		      bond_mat[i][i+NX-1]=1;
		      bond_mat[i][i+NX*2-1]=2;
		    }

		  
		  
		}
	      //ODD ROW
	      else
		{
		  if (i!=NX*(NY-1))
		    {
		      bond_mat[i][i+NX]=1;
		      bond_mat[i][i+1]=1;
		      bond_mat[i][i-NX]=1;
		      //PBC
		      bond_mat[i][i+NX-1]=1;
		      
		    }

		  else
		    {
		      //top left corner 
		      bond_mat[i][i+1]=1;
		      bond_mat[i][i-NX]=1;
		      bond_mat[i][i+NX-1]=1;
		      bond_mat[i][0]=1;
		      
		    }
		}
	    }
	  //RIGHT BOUNDARY  
	  else
	    {	
	      
		//EVEN ROW
	      //printf("\t%d", i);
	      if ((i/NX)%2==0)
		{
		  
		  //if it's not right corners 
		  if (i!=LEN-1 and i!=NX-1)
		    {
		    
		      bond_mat[i][i+NX]=1;
		      bond_mat[i][i-NX]=1;
		      bond_mat[i][i-1]=1;
		      //PBC
		      bond_mat[i][i-NX+1]=1;
		      
		    }
		  else
		    {
		      //printf("\n bottom right corner %d\n", i);
		      bond_mat[i][i+NX]=1;
		      bond_mat[i][i-1]=1;
		      //PBC
		      bond_mat[i][0]=1;
		      bond_mat[i][LEN-1]=1;
		      
		    }
		}
	      //ODD ROW
	      else
		{
		  //if not top right corner
		  if (i!=LEN-1)
		    {
		      
			bond_mat[i][i+NX]=1;
			bond_mat[i][i-NX]=1;
			bond_mat[i][i-NX-1]=2;
			bond_mat[i][i-1]=1;
			bond_mat[i][i+NX-1]=2;
			//PBC
			bond_mat[i][i+1]=2;
			bond_mat[i][i-NX+1]=1;
			bond_mat[i][i-NX-NX+1]=2;
			
		    }
		  else
		    {
		      //top right corner 
		      //printf("top right corner: %d\n",i);
		      bond_mat[i][i-NX]=1; 
		      bond_mat[i][i-NX-1]=2;
		      bond_mat[i][i-1]=1;
		      //PBC
		      bond_mat[i][0]=2;
		      bond_mat[i][i-NX+1]=1;
		      bond_mat[i][i-NX-NX+1]=2;
		      bond_mat[i][NX-1]=1;
		      bond_mat[i][NX-2]=2;
		    }
		}	
	      
	    }
	}
      //else if (i/NX==0)
      /* INTERMEDIATE + PBC IMAGES*/
      else
	{
	  /* ODD ROW and ODD COLUMN or EVEN ROW and EVEN COLUMN */
	  if (((i/NX)%2==0 and (i%NX)%2==0) or ((i/NX)%2==1 and (i%NX)%2==1))
	    {
	      //gurantee no problem as we exclude LEFT and RIGHT boundaries
	      bond_mat[i][i-1]=1;
	      bond_mat[i][i+1]=1;
	      
	      image_top = i+NX; 
	      if (image_top>LEN-1)
		image_top=i-NX*(NY-1);
	      image_bottom=i-NX;
	      if (image_bottom<0)
		image_bottom=i+NX*(NY-1);
	      
	      
	      bond_mat[i][image_top]=1;
	      bond_mat[i][image_top+1]=2;
	      bond_mat[i][image_top-1]=2;
	      bond_mat[i][image_bottom]=1;
	      bond_mat[i][image_bottom+1]=2;
	      bond_mat[i][image_bottom-1]=2;
	      

	      

	    }
	  /* EVEN ROW and ODD COLUMN or ODD ROW and EVEN COLUMN */
	  else
	    {
	      
	      //gurantee no problem as we exclude LEFT and RIGHT boundaries
	      bond_mat[i][i-1]=1;
	      bond_mat[i][i+1]=1;
	      image_top = i+NX;
	      if (image_top>LEN)
		image_top=i-NX*(NY-1);
	      
	      image_bottom=i-NX; 
	      if (image_bottom<0)
		image_bottom=i+NX*(NY-1);
	      
	      bond_mat[i][image_top]=1;
	      bond_mat[i][image_bottom]=1;
	      
	    }
	
	}
	
    }
  return 0;
}


/*	Bond Matrix is Symmetric	*/
/*	Function to check the same	*/



int check_bond_mat()
{
   for(int i=0;i<LEN;i++)
   {
	for(int j=0;j<i;j++)
	{
		if(bond_mat[i][j]!=bond_mat[j][i])
		{
			
			printf("At row = %d Col = %d\n",i,j);
			print_and_exit("ERROR: Bond Matrix is not Symmetric\n");
		}
		if(bond_mat[i][j] != 0)
		{
                        num_bonds++;
                }
	}
   }
   return 0;
}

/*	Printing the Bond Pairs, no repeats	*/
int bonds(FILE *fp)
{
   for(int i=0;i<LEN;i++)
   {
        for(int j=0;j<i;j++)
        {
		bondType=bond_mat[i][j]; 
		if(bondType!=0)
 	       	{   
			if (particle_id[i]==4 or particle_id[j]==4)
				{bondType=bondType+1;}
			else
				{bondType=bondType-1;}
			fprintf(fp,"%d,%d,%d\n",i,j,bondType);
	       		//FINAL output for bonds: 
			//0: normal horizontal/vertical 
	       		//1: normal diagonal 
	       		//2: dilated/shrinked horizontal/vertical 
	       		//3: dilated/shrinked diagonal  
		}
        }
   }
   return 0;
}

/*----------------------------------------------*/
/*	Generate dihedrals HOOMD convention	*/
/*----------------------------------------------*/
int generate_dihedrals()
{
   int j=0;
   int Nb_x = (NX-1)/2;
   int Nb_y = (NY-1)/2;
   int ii = 0;
   /*      Type I dihedral         */
   for (int nx=0;nx<Nb_x; nx++){
	for (int ny=0;ny<Nb_y; ny++)
	{
		
		ii =1+ NX*(ny*2+1)+nx*2;
		//set 'ii' as the center 4
		
		//TYPE A

		dihedrals[cnt_dihedrals][j]=ii-NX; j++;	//1
		dihedrals[cnt_dihedrals][j]=ii; j++;	//4
		dihedrals[cnt_dihedrals][j]=ii-NX-1; j++;//0
		dihedrals[cnt_dihedrals][j]=ii-1; j++;	//3
		cnt_dihedrals++;
		j=0;
		
		dihedrals[cnt_dihedrals][j]=ii+1; j++; //5
                dihedrals[cnt_dihedrals][j]=ii; j++;   //4
                dihedrals[cnt_dihedrals][j]=ii+NX+1; j++; //8
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;  //7
                cnt_dihedrals++;
		j=0;

		
 
		//TYPE B
		dihedrals[cnt_dihedrals][j]=ii-1; j++;  
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX-1; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;
                cnt_dihedrals++;
		j=0;


		dihedrals[cnt_dihedrals][j]=ii-NX; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii-NX+1; j++;
                dihedrals[cnt_dihedrals][j]=ii+1; j++;
                cnt_dihedrals++;
		j=0;

 
		//TYPE C
		dihedrals[cnt_dihedrals][j]=ii-NX-1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii-1; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX-1; j++;
                cnt_dihedrals++;
                j=0;
		
		dihedrals[cnt_dihedrals][j]=ii-NX-1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii-NX; j++;
                dihedrals[cnt_dihedrals][j]=ii-NX+1; j++;
                cnt_dihedrals++;
                j=0;

		dihedrals[cnt_dihedrals][j]=ii-NX+1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii+1; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX+1; j++;
                cnt_dihedrals++;
                j=0;
		
		dihedrals[cnt_dihedrals][j]=ii+NX+1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX-1; j++;
                cnt_dihedrals++;
                j=0;
	
		
	
	} }  
   for (int nx=0;nx<Nb_x-1; nx++){
        for (int ny=0;ny<Nb_y; ny++)
        {
		ii = 2+ NX*(ny*2+1)+nx*2;
		dihedrals[cnt_dihedrals][j]=ii-1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii-NX; j++;
                dihedrals[cnt_dihedrals][j]=ii+1; j++;
                cnt_dihedrals++;
                j=0;

		dihedrals[cnt_dihedrals][j]=ii-1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;
                dihedrals[cnt_dihedrals][j]=ii+1; j++;
                cnt_dihedrals++;
                j=0;

		

	}}
   for (int nx=0;nx<Nb_x; nx++){
        for (int ny=0;ny<Nb_y-1; ny++)
        {	
		ii =  1+ NX*(ny*2+2)+nx*2;
		dihedrals[cnt_dihedrals][j]=ii-NX; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii-1; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;
                cnt_dihedrals++;
                j=0;

                dihedrals[cnt_dihedrals][j]=ii-NX; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii+1; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;
                cnt_dihedrals++;
                j=0;

	}}
   printf ("#dihedrals: %d", cnt_dihedrals);
   return 0;
}	



int generate_dihedrals_test()
{
   int j=0;
   int Nb_x = (NX-1)/2;
   int Nb_y = (NY-1)/2;
   int ii = 0;
   /*      Type I dihedral         */
   for (int nx=0;nx<Nb_x; nx++){
	for (int ny=0;ny<Nb_y; ny++)
	{
		
		ii =1+ NX*(ny*2+1)+nx*2;
		//set 'ii' as the center 4
		
		//TYPE A

		dihedrals[cnt_dihedrals][j]=ii-NX; j++;	//1
		dihedrals[cnt_dihedrals][j]=ii; j++;	//4
		dihedrals[cnt_dihedrals][j]=ii-NX-1; j++;//0
		dihedrals[cnt_dihedrals][j]=ii-1; j++;	//3
		cnt_dihedrals++;
		j=0;
		
		dihedrals[cnt_dihedrals][j]=ii+1; j++; //5
                dihedrals[cnt_dihedrals][j]=ii; j++;   //4
                dihedrals[cnt_dihedrals][j]=ii+NX+1; j++; //8
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;  //7
                cnt_dihedrals++;
		j=0;

		
 
		//TYPE B
		dihedrals[cnt_dihedrals][j]=ii-1; j++;  
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX-1; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;
                cnt_dihedrals++;
		j=0;


		dihedrals[cnt_dihedrals][j]=ii-NX; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii-NX+1; j++;
                dihedrals[cnt_dihedrals][j]=ii+1; j++;
                cnt_dihedrals++;
		j=0;

 
		//TYPE C
		dihedrals[cnt_dihedrals][j]=ii-NX-1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii-1; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX-1; j++;
                cnt_dihedrals++;
                j=0;
		
		dihedrals[cnt_dihedrals][j]=ii-NX-1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii-NX; j++;
                dihedrals[cnt_dihedrals][j]=ii-NX+1; j++;
                cnt_dihedrals++;
                j=0;

		dihedrals[cnt_dihedrals][j]=ii-NX+1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii+1; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX+1; j++;
                cnt_dihedrals++;
                j=0;
		
		dihedrals[cnt_dihedrals][j]=ii+NX+1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX-1; j++;
                cnt_dihedrals++;
                j=0;
	
		
	
	} }  
   for (int nx=0;nx<Nb_x-1; nx++){
        for (int ny=0;ny<Nb_y; ny++)
        {
		ii = 2+ NX*(ny*2+1)+nx*2;
		dihedrals[cnt_dihedrals][j]=ii-1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii-NX; j++;
                dihedrals[cnt_dihedrals][j]=ii+1; j++;
                cnt_dihedrals++;
                j=0;

		dihedrals[cnt_dihedrals][j]=ii-1; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;
                dihedrals[cnt_dihedrals][j]=ii+1; j++;
                cnt_dihedrals++;
                j=0;

		

	}}
   for (int nx=0;nx<Nb_x; nx++){
        for (int ny=0;ny<Nb_y-1; ny++)
        {	
		ii =  1+ NX*(ny*2+2)+nx*2;
		dihedrals[cnt_dihedrals][j]=ii-NX; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii-1; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;
                cnt_dihedrals++;
                j=0;

                dihedrals[cnt_dihedrals][j]=ii-NX; j++;
                dihedrals[cnt_dihedrals][j]=ii; j++;
                dihedrals[cnt_dihedrals][j]=ii+1; j++;
                dihedrals[cnt_dihedrals][j]=ii+NX; j++;
                cnt_dihedrals++;
                j=0;

	}}
   printf ("#dihedrals: %d", cnt_dihedrals);
   return 0;
}	



int generate_dihedrals_pbc()
{
   int j=0;
   int Nb_x = (NX-2)/2;
   int Nb_y = (NY-2)/2;
   int ii = 0;
   
   int i0=0, i1=0, i2=0, i3=0, i4=0, i5=0, i6=0, i7=0, i8=0;
   
    /* 
       x: left/right boundary 
       y: top/bottom boundary
       *: intermediate 
       
       TYPE I 
       y y y y y y c
       * * * * * * x
       * * * * * * x
       * * * * * * x 
       * * * * * * x
       * * * * * * x
       
       TYPE II 

       TYPE III
       
       
    */
   
   
   
   /*      Type I dihedral INTERMEDIATE         */
   for (int nx=0;nx<=Nb_x; nx++){
     for (int ny=0;ny<=Nb_y; ny++)
       {
		
	 ii =1+ NX*(ny*2+1)+nx*2;
	 //set 'ii' as the center 4
	 
	 //left boundary to intermediate
	 if (nx<Nb_x and ny<Nb_y)
	   {
	     i4 = ii; 
	     i0 = i4-NX-1;
	     i1 = i4-NX;
	     i2 = i4-NX+1;
	     i3 = i4-1;
	     i5 = i4+1;
	     i6 = i4+NX-1; 
	     i7 = i4+NX; 
	     i8 = i4+NX+1;
	   }
	 //right boundary 
	 else if (nx==Nb_x and ny<Nb_y)
	   {
	     i4 = ii; 
	     i0 = i4-NX-1;
	     i1 = i4-NX;
	     i2 = i4+1-NX-NX; //pbc
	     i3 = i4-1;
	     i5 = i4+1-NX;//pbc
	     i6 = i4+NX-1; 
	     i7 = i4+NX; 
	     i8 = i4+1; //pbc
	   }
	 //top boundary 
	 else if (nx<Nb_x and ny==Nb_y)
	   {
	     i4 = ii; 
	     i0 = i4-NX-1;
	     i1 = i4-NX;
	     i2 = i4-NX+1;
	     i3 = i4-1;
	     i5 = i4+1;
	     i6 = i4+NX-1-NX*NY; //pbc 
	     i7 = i4+NX-NX*NY; //pbc
	     i8 = i4+NX+1-NX*NY;//pbc
	   }
	 //top right corner 
	 else
	   {
	     i4 = ii; 
	     i0 = i4-NX-1;
	     i1 = i4-NX;
	     i2 = i4-NX+1-NX; //PBC
	     i3 = i4-1;
	     i5 = i4+1-NX; //PBC
	     i6 = i4+NX-1-NX*NY; //PBC 
	     i7 = i4+NX-NX*NY; //PBC
	     i8 = 0;//PBC
      
	   }

	 
	 //TYPE A
	 
	 dihedrals[cnt_dihedrals][j]=i1; j++;
	 dihedrals[cnt_dihedrals][j]=i4; j++;
	 dihedrals[cnt_dihedrals][j]=i0; j++;
	 dihedrals[cnt_dihedrals][j]=i3; j++;
	 cnt_dihedrals++;
	 j=0;
	 
	 dihedrals[cnt_dihedrals][j]=i5; j++;
	 dihedrals[cnt_dihedrals][j]=i4; j++;
	 dihedrals[cnt_dihedrals][j]=i8; j++;
	 dihedrals[cnt_dihedrals][j]=i7; j++;
	 cnt_dihedrals++;
	 j=0;
	 
	 
 
	 //TYPE B
	 dihedrals[cnt_dihedrals][j]=i3; j++;
	 dihedrals[cnt_dihedrals][j]=i4; j++;
	 dihedrals[cnt_dihedrals][j]=i6; j++;
	 dihedrals[cnt_dihedrals][j]=i7; j++;
	 cnt_dihedrals++;
	 j=0;
	 
	 
	 dihedrals[cnt_dihedrals][j]=i1; j++;
	 dihedrals[cnt_dihedrals][j]=i4; j++;
	 dihedrals[cnt_dihedrals][j]=i2; j++;
	 dihedrals[cnt_dihedrals][j]=i5; j++;
	 cnt_dihedrals++;
	 j=0;
	 
 
	 //TYPE C
	 dihedrals[cnt_dihedrals][j]=i0; j++;
	 dihedrals[cnt_dihedrals][j]=i4; j++;
	 dihedrals[cnt_dihedrals][j]=i3; j++;
	 dihedrals[cnt_dihedrals][j]=i6; j++;
	 cnt_dihedrals++;
	 j=0;
	 
	 dihedrals[cnt_dihedrals][j]=i0; j++;
	 dihedrals[cnt_dihedrals][j]=i4; j++;
	 dihedrals[cnt_dihedrals][j]=i1; j++;
	 dihedrals[cnt_dihedrals][j]=i2; j++;
	 cnt_dihedrals++;
	 j=0;
	 
	 dihedrals[cnt_dihedrals][j]=i2; j++;
	 dihedrals[cnt_dihedrals][j]=i4; j++;
	 dihedrals[cnt_dihedrals][j]=i5; j++;
	 dihedrals[cnt_dihedrals][j]=i8; j++;
	 cnt_dihedrals++;
	 j=0;
	 
	 dihedrals[cnt_dihedrals][j]=i8; j++;
	 dihedrals[cnt_dihedrals][j]=i4; j++;
	 dihedrals[cnt_dihedrals][j]=i7; j++;
	 dihedrals[cnt_dihedrals][j]=i6; j++;
	 cnt_dihedrals++;
	 j=0;
	
		
	
       } 
   }  



   


   /* Type II dihedral neighboring along x */
   for (int nx=0;nx<=Nb_x; nx++){
     for (int ny=0;ny<=Nb_y; ny++)
       {
	 
	 
	 
	 //intermediate
	 if (nx<Nb_x and ny<Nb_y)
	   {
	     ii = 2+ NX*(ny*2+1)+nx*2;
	     i4=ii;
	     i0 = i4-NX-1;
	     i1 = i4-NX;
	     i2 = i4-NX+1;
	     i3 = i4-1;
	     i5 = i4+1;
	     i6 = i4+NX-1; 
	     i7 = i4+NX; 
	     i8 = i4+NX+1;
	     
	     
	   }
	 
	 //right boundary connects to left boundary excluding corner 
	 else if (nx==Nb_x and ny<Nb_y)
	   {
	     ii = 2+ NX*(ny*2+1)+nx*2-NX;//shift to the right by NX look for PBC (right boundry)
	     i4 = ii;
	     i0 = i4-NX-1+NX;//pbc
	     i1 = i4-NX;
	     i2 = i4-NX+1-NX;
	     i3 = i4-1+NX;//pbc
	     i5 = i4+1;
	     i6 = i4+NX-1+NX; //pbc
	     i7 = i4+NX; 
	     i8 = i4+NX+1;
	   }
	 
	 //top boundary excluding top right corner
	 else if (nx<Nb_x and ny==Nb_y)
	   {
	     
	     ii = 2+ NX*(ny*2+1)+nx*2;
	     i4 = ii;
	     i0 = i4-NX-1;
	     i1 = i4-NX;
	     i2 = i4-NX+1;
	     i3 = i4-1;
	     i5 = i4+1;
	     i6 = i4+NX-1-NX*NY; 
	     i7 = i4+NX-NX*NY; 
	     i8 = i4+NX+1-NX*NY;
	   }
	 
	 //top right corner (NOT DONE)
	 else
	   {
	     
	     ii = 2+ NX*(ny*2+1)+nx*2-NX; //shifted to the left 
	     i4 = ii;
	     i0 = i4-NX-1+NX;//pbc
	     i1 = i4-NX;
	     i2 = i4-NX+1;
	     i3 = i4-1;
	     i5 = i4+1;
	     i6 = 0+NX-1; //pbc corner
	     i7 = 0; //pbc
	     i8 = 1;//pbc
	   }

	 dihedrals[cnt_dihedrals][j]=i3; j++;//3
	 dihedrals[cnt_dihedrals][j]=i4; j++;//4
	 dihedrals[cnt_dihedrals][j]=i1; j++;//1
	 dihedrals[cnt_dihedrals][j]=i5; j++;//5
	 cnt_dihedrals++;
	 j=0;
	 
	 dihedrals[cnt_dihedrals][j]=i3; j++;//3
	 dihedrals[cnt_dihedrals][j]=i4; j++;//4
	 dihedrals[cnt_dihedrals][j]=i7; j++;//7
	 dihedrals[cnt_dihedrals][j]=i5; j++;//5
	 cnt_dihedrals++;
	 j=0;
	 
	 
	 
       }
   }
   for (int nx=0;nx<=Nb_x; nx++){
     for (int ny=0;ny<=Nb_y; ny++)
       {	
	 
	 
	  //intermediate
	 if (nx<Nb_x and ny<Nb_y)
	   {
	     ii = 1+ NX*(ny*2+2)+nx*2;
	     i4=ii;
	     i0 = i4-NX-1;
	     i1 = i4-NX;
	     i2 = i4-NX+1;
	     i3 = i4-1;
	     i5 = i4+1;
	     i6 = i4+NX-1; 
	     i7 = i4+NX; 
	     i8 = i4+NX+1;
	     
	     
	   }
	  //right boundary connects to left boundary excluding corner 
	 else if (nx==Nb_x and ny<Nb_y)
	   {
	     ii = 1+ NX*(ny*2+2)+nx*2;
	     i4=ii;
	     i0 = i4-NX-1;
	     i1 = i4-NX;
	     i2 = i4-NX+1-NX;//pbc
	     i3 = i4-1;
	     i5 = i4+1-NX;//pbc
	     i6 = i4+NX-1; 
	     i7 = i4+NX; 
	     i8 = i4+NX+1-NX;
	     
	   }
	 
	 else if (nx<Nb_x and ny==Nb_y)
	   {
	     ii = 1+ NX*(ny*2+2)+nx*2-NX*NY; //shifted to the bottom
	     i4=ii;
	     i0 = i4-NX-1 + NX*NY;//pbc
	     i1 = i4-NX + NX*NY;//pbc
	     i2 = i4-NX+1 + NX*NY;//pbc
	     i3 = i4-1;
	     i5 = i4+1;
	     i6 = i4+NX-1; 
	     i7 = i4+NX; 
	     i8 = i4+NX+1;
	     
	   }
	 //top right corner 
	 else
	   {
	   
	     ii = 1+ NX*(ny*2+2)+nx*2-NX*NY; //shifted to the bottom
	     
	     i4=ii;
	     i0 = i4-NX-1 + NX*NY;//pbc
	     i1 = i4-NX + NX*NY;//pbc
	     i2 = i4-NX+1 + NX*NY-NX;//pbc
	     i3 = i4-1;
	     i5 = i4+1-NX;//pbc
	     i6 = i4+NX-1; 
	     i7 = i4+NX; 
	     i8 = i4+NX+1-NX;//pbc

	     
	   }

	 dihedrals[cnt_dihedrals][j]=i1; j++;//1
	 dihedrals[cnt_dihedrals][j]=i4; j++;//4
	 dihedrals[cnt_dihedrals][j]=i3; j++;//3
	 dihedrals[cnt_dihedrals][j]=i7; j++;//7
	 cnt_dihedrals++;
	 j=0;
	 
	 dihedrals[cnt_dihedrals][j]=i1; j++;//1
	 dihedrals[cnt_dihedrals][j]=i4; j++;//4
	 dihedrals[cnt_dihedrals][j]=i5; j++;//5
	 dihedrals[cnt_dihedrals][j]=i7; j++;//7
	 cnt_dihedrals++;
	 j=0;
	 
       }
   }
   printf ("#dihedrals: %d", cnt_dihedrals);
   return 0;
}	



/*		Print dihedrals		*/		
int out_dihedrals(FILE *fp)
{
   for(int i=0;i<cnt_dihedrals;i++)
   {
	fprintf(fp,"%d,%d,%d,%d\n",dihedrals[i][0],dihedrals[i][1],dihedrals[i][2],dihedrals[i][3]);
   }
   return 0;
}

/*	Particle Typeid		*/
/*	Normal nodes - particle Id 0 */
/*	Clamped left - particle Id  1	*/
/*	Right end free to slide along X - - particle Id  3 */
/*	Backbone of ribbon - particle Id 4 */

int particle_typeid()
{
   for(int i=0;i<LEN;i++)
   {
	if(i%NX==0 || i%NX==1) //Clamping two columns of lattice sites on the left
	{
		particle_id[i]=1;
		//printf("particle_id[%d] = %d\n",i,particle_id[i]);
	}
	else if (i%NX==NX-1 || i%NX==NX-2) //Two cols lattice sites on the right constrained to move only X
	{
		particle_id[i]=1;
		//printf("particle_id[%d] = %d\n",i,particle_id[i]);
	}
	else
                particle_id[i]=0;//Normal lattice sites 
  }
  //Backbone of the ribbon excluding two lattice sites at each boundary
  if(defect==1)  //line defect along the backbone 
	for(int i=2;i<NX-2;i++)
  	{	if (i%spx==0)
		{	particle_id[(NY/2)*NX + i]=4; //bumps
		}
		else	
			particle_id[(NY/2)*NX + i]=0; //normal
  	}
  else if (defect==2)
    for(int i=edgex; i<NX-edgex; i+=spx)
      {	
	for (int j=edgey; j<NY-edgey; j+=spy)
	  
	  {       
	    particle_id[j*NX + i]=4; 
	    //ip->z[j*NX+i]=0.1*pow(-1, i/4+j/4);	//bumps
	  }
	
      }

else if (defect==-2) //for negative omega, still same particle type as defect=2
    for(int i=edgex; i<NX-edgex; i+=spx)
      {
        for (int j=edgey; j<NY-edgey; j+=spy)

          {       
            particle_id[j*NX + i]=4; 
            
	}	
	}
	//ip->z[j*NX+i]=0.1*pow(-1, i/4+j/4);       //bumps
            //          }
            //
            //                }
            //


   //else if (defect==2)
   //     for(int i=3; i<NX-3;i++)
   //     {       for (int j=3; j<NY-3; j++)
   //                    if (i%spx==0 and j%spy==0)
   //                     {       particle_id[j*NX + i]=4;
                                //ip->z[j*NX+i]=0.1*pow(-1, i/4+j/4);   //bumps
                                //                        }
                                //
                                //                                }
                                //
  //(n,n) tiling 
  else if (defect==3)
    for(int i=2; i<NX-2;i++)
      {       for (int j=2; j<NY; j+=2)
	  if(j%(spy*2)==0)
	    {
	      if (i%(spx*3)==0)
		{       particle_id[j*NX + i]=4; //bumps
		}
	    }
	  else
	    {
	      if ((i-3)%(spx*3)==0)
		{       particle_id[j*NX + i]=4; //bumps
		}
	      
	    }

      }
  else if (defect==4)
	{
	  for (int i=0; i<LEN; i++)
	    {
	      if (i==NX*(NY/2)+NX/2)
		{particle_id[i]=4;}
	    }
	}
  

  
   return 0;
}	

/*	Print particles TypeId		*/
int out_typeId(FILE *fp)
{
   for(int i=0;i<LEN;i++)
   {
	fprintf(fp,"%u\n",particle_id[i]);
   }
   return 0;
}

void print_and_exit(char *format, ...)
{
    va_list list;
    va_start(list,format);
    vfprintf(stderr,format,list);
    va_end(list);
    exit(1);
}

