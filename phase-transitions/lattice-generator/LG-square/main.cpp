/*
 * Code to generate sqaure lattice with defects 
 * Author: Paul Hanakata, Harvard University 2020-2022
 * For more details go to: https://github.com/phanakata/programmable-matter
 * See our publications: 
 * PZ Hanakata, A Plummer, DR Nelson PRL 128 (7), 075902
 * A Plummer, PZ Hanakata, DR Nelson PRM 6 (11), 115203
 * */
#include <stdio.h>
#include <stdlib.h>
//#include <cstdlib>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <limits>
//#include <iostream>
//#include <fstream>
//#include <string>
//#include <sstream>
//#include <iomanip>
//#include <vector>
#include "lattice_variables.h"
#include "lattice.h"
#include "streamCompaction.h"

int NX,NY,LEN,defect,spx, spy, pbc, edgex, edgey;

int main( int argc, char **argv )
{

   switch (argc){
     case 7:
       sscanf(argv[1],"%d",&NX); //number of nodes in x-direction
       sscanf(argv[2],"%d",&NY); //number of nodes in y-direction
       sscanf(argv[3],"%d",&pbc); //pbc:0 (FALSE), pbc:1 (TRUE)
       sscanf(argv[4],"%d",&defect); //defect:2/-2 postive/negative defect in 2d arrays, defect:4 defect at the center 
       sscanf(argv[5],"%d",&spx); //spacing between defects in x-direction
       sscanf(argv[6],"%d",&spy); //spacing between defects in y-direciton 
       break;
     default:
       print_and_exit("Usage: %s NX NY pbc defect_type spacing_x spacing_y\n",argv[0]);
   }

   LEN = NX*NY;
   if (pbc==0)
	{
	edgex = 6;
	edgey = 6; 
	}
   else
	{
	edgex = spx-spx;  
        edgey = spy-spy;
	}
   /*	Output File	*/
   FILE *lat;

   /*      Character array for directory pathname and filename     */
   char filepath[256];

   /*      Filepaths for output files      */
   //printf("Enter path where the lattice.dat file will be written\n");
   sprintf(filepath,"../configurations/lattice_pbc%d_d%d_L%d_W%d_spx%d_spy%d.dat", pbc, defect,NX,NY,spx,spy);
   printf("Filename of Lattice Details: %s\n",filepath);
   lat = fopen(filepath, "w");
   if (lat == NULL)
   {
   	print_and_exit("Could Not Open File:lattice.dat");
   }


   /*      Byte size of structures         */
   size_t nBytes = sizeof(latticeStruct); 
   
   /*      Dynamic memory allocation for coordinates       */
   latticeStruct *h_coords = (latticeStruct *)malloc(nBytes);

   /*      Lattice positions	    */
   initialLatticeStruct(h_coords, LEN);
  
   /*      Generating the Bond Matrix	*/
   if (pbc==0)
     {
       lattice_connectivity();
     }
   else
     {
       lattice_connectivity_pbc();
     }

/*   for(int i=0;i<LEN;i++)
   {
	for(int j=0; j<LEN;j++)
	{
		printf("%d ",bond_mat[i][j]);
	}
	printf("\n");
   }
*/   
   /*	Checking if Bond Matrix is Symmetric	*/
   check_bond_mat(); 

   /*	Generating the Dihedrals	*/
   if (pbc==0)
     {
       generate_dihedrals();
     }
   else
     {
       generate_dihedrals_pbc();

     }
   /*	Sorting Dihedrals using 2nd particle as pivot	*/
   //insertionSortDihedrals(cnt_dihedrals);

   /*	Particle Type ID	*/
   particle_typeid();

   /*	Printing lattice configuration	*/
   /*	Total Particles		*/
   fprintf(lat,"%d\n",LEN);
   printf("#particles %d\n",LEN);

   /*	Initial Paerticle Position in Flat configuration	*/
   for(int i=0;i<LEN;i++)
   {
        fprintf(lat,"%.8f,%.8f,%.8f\n",h_coords->x[i],h_coords->y[i],h_coords->z[i]);
   }
   /*   Printing the Bond pairs */
   fprintf(lat,"%d\n",num_bonds);
   printf("#Bonds %d\n",num_bonds);
   bonds(lat);

   /*   Printing the Dihedrals  */
   fprintf(lat,"%d\n",cnt_dihedrals);
   printf("#dihedrals %d\n",cnt_dihedrals);
   out_dihedrals(lat);

   /*   Printing particle type Ids	*/
   out_typeId(lat);

   /*	Close output file	*/
   fclose(lat);

   bond_compaction();
   return 0;
}

