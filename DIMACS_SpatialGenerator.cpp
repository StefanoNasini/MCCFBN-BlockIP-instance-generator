
/**************************************************************************************
 * DIMACS Minimum Cost Flow in Bipartite Graphs 
 **************************************************************************************/

#include <iostream>
#include <fstream>
using namespace std;
#include <cstring>
#include <cctype>
#include "math.h"
#include <vector>
#include <cstdlib>

int a = 1664525;
int c = 1013904223;
int mm = 4294967295;

unsigned int My_Rand()
{
    // My initial seed is 21051983
    static unsigned int seed = 21051983;
 
    // Take the current seed and generate a new value from it
    seed = (a * seed + c);
 
    // Take the seed and return a value between 0 and m
    return seed  % mm;
}







int main(int argc, char* argv[]) {

	if (argc < 7 ) {
		cerr << "usage: DIMACS_Generator out_filename ORIG DEST ALPHA SpaceDist MM FRACTIONAL [factor]\n";
		exit(-1);
	}

	// -------------------------------------------------------------------------------
	// INPUT PARAMETERS
	// -------------------------------------------------------------------------------
	// - ORIG and DEST are the primary and secondary layers of the bipartite graph. 
	// - The parameter AA (in [0,1]) controls the angle between the two lines of origens
	//	and destinations.
	// - The parameters SpaceDist is used to select the spatial distribution of points.
	// - The parameters MM (in [0,1]) controls the size of random shocks of points around
	//	their deterministic coordinates. 
	// - The parameters FRACTIONAL is an integer number which denotes whether the problem 
	//	parameters are supposed to be integer or allowed to be fractional. 
	// - The optional parameter factor (in [0,1]) is used to produce excess capacity. The closer
	//	to 0 the larger is the excess capacity.
	// -------------------------------------------------------------------------------
	
	int ORIG;		// number of facilities;
	int DEST;		// number of destination;
	double AA;		// angle in radiant between the two lines of facilities and destinations (in [0,1]);
	int SpaceDist;		// type of spatial distribution (either linear or spherical)
	double MM;
	double eps = 0.001;	
	double factor; //= 0.99;
	int FRACTIONAL;

	sscanf(argv[2], "%i", &ORIG);
	sscanf(argv[3], "%i", &DEST);
	sscanf(argv[4], "%lf",&AA);
	sscanf(argv[5], "%i", &SpaceDist);
	sscanf(argv[6], "%lf", &MM);
	sscanf(argv[7], "%i", &FRACTIONAL);
	if (argc >= 9 ) {
		sscanf(argv[8], "%lf", &factor );
	}else{
		factor = (double)0.99;
	}

	ofstream MyOutputFile;
	MyOutputFile.open (argv[1]);

	// DIMACS: http://lpsolve.sourceforge.net/5.5/DIMACS_maxf.htm

	MyOutputFile << "c This is the data input of a minimum cost flow problem in bipartite graphs in DIMACS format." << endl;
	MyOutputFile << "c The problem declaration consist in specifying the the number of nodes in the primary and secondary layers and the number of arcs." << endl;
	MyOutputFile << "c format: p bip NODES_1 NODES_2 ARCS" << endl;
	MyOutputFile << "c format: n1 ID WHICH" << endl;
	MyOutputFile << "c format: n2 ID WHICH" << endl;
	MyOutputFile << "c format: a SRC DST LOWER-CAP UPPER-CAP COST" << endl;
	MyOutputFile << "c ----------------------------------" << endl;

	double CostShortage = 100000;
	double Scal = DEST + ORIG;
	int Strech = round((double)DEST/((double)ORIG));

	try {

		MyOutputFile << "c ---------------------------------- Problem line (nodes, links)" << endl;
		MyOutputFile << "p bip " << ORIG << " " << DEST << " " << ORIG*DEST << endl;

		if(FRACTIONAL < 1 ){ //--------------------------------------------------------------------

			//---------------------------------------------------------------------------------
			///////////////// INITIALIZATION //////////////////////////////////////////////////
			//---------------------------------------------------------------------------------

			int *supply = new int[ORIG];
			int *demand = new int[DEST];
			int *c = new int [DEST*ORIG];	// linear costs origin-destination
			int *ub = new int[DEST*ORIG];	// upper bounds on the decision variables
			double total_supply;
			int Large = ORIG*DEST*10;

			int k;
			double ran_x;
			double ran_y;

			//---------------------------------------------------------------------------------
			///////////////// DATA FILE BUILDING /////////////////////////////////////////////
			//---------------------------------------------------------------------------------

			if(SpaceDist == 1){

			//---------------------------------------------------------------------------------
			// Linear spatial distribution with linearly growing supplies and uniform demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(ORIG*(1+i) + max(0,DEST/(ORIG) - ORIG*(ORIG+1)/2));
				  total_supply = total_supply + supply[i];
				}

				cout << 0.1*(1-factor)*total_supply/ORIG << "\n";

				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(supply[i] + 0.1*(1-factor)*total_supply/ORIG);
				  MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
				  demand[j] = (int)max((int)round(factor*(total_supply- eps)/DEST),1);
				  //demand[j] = max((int)round((total_supply/DEST)),1);
				  MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        		c[k] = (abs(Strech*(i+1)-((int)round(ran_x*j)+1)) + abs(Strech*(i+1) - abs(((int)round(ran_y*j)+1)-(int)round(AA*DEST))));
						//ub[k] = (int)round(demand[j]*(1/(ORIG) + c[k]/(Scal)));
						ub[k] = (int)round(demand[j]*(1 + c[k]/(Scal)));
						c[k] = (int)round(0.5 + c[k] % 10 + 2*min(DEST,ORIG)*c[k]/max(DEST,ORIG) );

					  	MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
					  	k++;

					}
				}

			} else if(SpaceDist == 2){

			//---------------------------------------------------------------------------------
			// Linear spatial distribution with quadratically growing supplies and uniform demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(ORIG*(1+i)*(1+i) + max(0,DEST/ORIG - ORIG*ORIG*ORIG/3 - ORIG*ORIG/2 - ORIG/6));
				  total_supply = total_supply + supply[i];
				}

				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(supply[i] + 0.1*(1-factor)*total_supply/ORIG);
				  MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
				  demand[j] = (int)max((int)round(factor*(total_supply- eps)/DEST),1);
				  //demand[j] = max((int)round((factor*total_supply/DEST)),1);
				  MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-((int)round(ran_x*j)+1)) + abs(Strech*(i+1) - abs(((int)round(ran_y*j)+1)-(int)round(AA*DEST))));
						//ub[k] = (int)round(demand[j]*(1/(ORIG) + c[k]/(Scal)));
						ub[k] = (int)round(demand[j]*(1 + c[k]/(Scal)));

						c[k] = (int)round(0.5 + c[k] % 10 + 2*min(DEST,ORIG)*c[k]/max(DEST,ORIG) );

					 	MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
					 	k++;

					}
				}

			} else if(SpaceDist == 3){

			//---------------------------------------------------------------------------------
			// Circular distribution of customers with linearly growing supplies and uniform demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(ORIG*(1+i) + max(0,DEST/ORIG - ORIG*(ORIG+1)/2));
				  total_supply = total_supply + supply[i];
				}

				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(supply[i] + 0.1*(1-factor)*total_supply/ORIG);
				  MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
				  demand[j] = (int)max((int)round(factor*(total_supply- eps)/DEST),1);
				  //demand[j] = max((int)round((factor*total_supply/DEST)),1);
				  MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-(int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_x*j)/DEST)/2)) + abs(Strech*(i+1) - (int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_y*j)/DEST)/2)) );
						//ub[k] = (int)round(demand[j]*(1/(ORIG) + c[k]/(Scal)));
						ub[k] = (int)round(demand[j]*(1 + c[k]/(Scal)));

						c[k] = (int)round(0.5 + c[k] % 10 + 2*min(DEST,ORIG)*c[k]/max(DEST,ORIG) );

					  	MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
					  	k++;
					  
					}
				}

			} else if(SpaceDist == 4){

			//---------------------------------------------------------------------------------
			// Circular distribution of customers with quadratically growing supplies and uniform demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(ORIG*(1+i)*(1+i) + max(0,DEST/ORIG - ORIG*ORIG*ORIG/3 - ORIG*ORIG/2 - ORIG/6));
				  total_supply = total_supply + supply[i];
				}

				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(supply[i] + 0.1*(1-factor)*total_supply/ORIG);
				  MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
				  demand[j] = (int)max((int)round(factor*(total_supply- eps)/DEST),1);
				  //demand[j] = max((int)round((factor*total_supply/DEST)),1);
				  MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        		c[k] = (abs(Strech*(i+1)-(int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_x*j)/DEST)/2)) + abs(Strech*(i+1) - (int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_y*j)/DEST)/2)) );
						//ub[k] = (int)round(demand[j]*(1/(ORIG) + c[k]/(Scal)));
						ub[k] = (int)round(demand[j]*(1 + c[k]/(Scal)));

						c[k] = (int)round(0.5 + c[k] % 10 + 2*min(DEST,ORIG)*c[k]/max(DEST,ORIG) );
				  
					  	MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
					  	k++;

					}
				}

			} else if(SpaceDist == 5){

			//---------------------------------------------------------------------------------
			// Linear spatial distribution with linearly growing supplies and linearly growing demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(ORIG*(1+i) + max(0,DEST/ORIG - ORIG*(ORIG+1)/2));
				  total_supply = total_supply + supply[i];
				}

				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(supply[i] + 0.1*(1-factor)*total_supply/ORIG);
				  MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
				  demand[DEST-j-1] = max((int)round(((factor*2*total_supply*(j + 1))/(DEST*(DEST+1)))),1);
				  MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-((int)round(ran_x*j)+1)) + abs(Strech*(i+1) - abs(((int)round(ran_y*j)+1)-(int)round(AA*DEST))));
						//ub[k] = (int)round(demand[j]*(1/(ORIG) + c[k]/(Scal)));
						ub[k] = (int)round(demand[j]*(1 + c[k]/(Scal)));

						c[k] = (int)round(0.5 + c[k] % 10 + 2*min(DEST,ORIG)*c[k]/max(DEST,ORIG) );

					  	MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
					 	k++;

					}
				}

			} else if(SpaceDist == 6){

			//---------------------------------------------------------------------------------
			// Linear spatial distribution with quadratically growing supplies and linearly growing demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(ORIG*(1+i)*(1+i) + max(0,DEST/ORIG - ORIG*ORIG*ORIG/3 - ORIG*ORIG/2 - ORIG/6));
				  total_supply = total_supply + supply[i];
				}

				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(supply[i] + 0.1*(1-factor)*total_supply/ORIG);
				  MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
				  demand[j] = max((int)round(((factor*2*total_supply*(j + 1))/(DEST*(DEST+1)))),1);
				  MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] =  (abs(Strech*(i+1)-((int)round(ran_x*j)+1)) + abs(Strech*(i+1) - abs(((int)round(ran_y*j)+1)-(int)round(AA*DEST))));
						//ub[k] = 1 + (int)round(demand[j]*(1/(ORIG) + c[k]/(Scal)));
						ub[k] = (int)round(demand[j]*(1 + c[k]/(Scal)));

						c[k] = (int)round(0.5 + c[k] % 10 + 2*min(DEST,ORIG)*c[k]/max(DEST,ORIG) );

					 	MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
					 	k++;

					}
				}

			} else if(SpaceDist == 7){

			//---------------------------------------------------------------------------------
			// Circular distribution of customers with linearly growing supplies and linearly growing demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(ORIG*(1+i) + max(0,DEST/ORIG - ORIG*(ORIG+1)/2));
				  total_supply = total_supply + supply[i];
				}

				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(supply[i] + 0.1*(1-factor)*total_supply/ORIG);
				  MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
				  demand[j] = max((int)round(((factor*2*total_supply*(j + 1))/(DEST*(DEST+1)))),1);
				  MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        		c[k] = (abs(Strech*(i+1)-(int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_x*j)/DEST)/2)) + abs(Strech*(i+1) - (int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_y*j)/DEST)/2)) );
						//ub[k] = 1 + (int)round(demand[j]*(1/(ORIG) + c[k]/(Scal)));
						ub[k] = (int)round(demand[j]*(1 + c[k]/(Scal)));

						c[k] = (int)round(0.5 + c[k] % 10 + 2*min(DEST,ORIG)*c[k]/max(DEST,ORIG) );

					  	MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
					  	k++;

					}
				}

			} else if(SpaceDist == 8){

			//---------------------------------------------------------------------------------
			// Circular distribution of customers with quadratically growing supplies and linearly growing demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(ORIG*(1+i)*(1+i) + max(0,DEST/ORIG - ORIG*ORIG*ORIG/3 - ORIG*ORIG/2 - ORIG/6));
				  total_supply = total_supply + supply[i];
				}

				for(int i=0; i<ORIG; i++){
				  supply[i] = (int)round(supply[i] + 0.1*(1-factor)*total_supply/ORIG);
				  MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
				  demand[j] = max((int)round(((factor*2*total_supply*(j + 1))/(DEST*(DEST+1)))),1);
				  MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-(int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_x*j)/DEST)/2)) + abs(Strech*(i+1) - (int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_y*j)/DEST)/2)) );
						//ub[k] = 1 + (int)round(demand[j]*(1/(ORIG) + c[k]/(Scal)));
						ub[k] = (int)round(demand[j]*(1 + c[k]/(Scal)));

						c[k] = (int)round( 0.5 + c[k] % 10 + 2*min(DEST,ORIG)*c[k]/max(DEST,ORIG) );

					  	MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
					  	k++;

					}
				}
			}

			delete[] c;
			delete[] ub;
			delete[] supply;
			delete[] demand;  

		}else{ //----------------------------------------------------------------------------------

			//---------------------------------------------------------------------------------
			///////////////// INITIALIZATION //////////////////////////////////////////////////
			//---------------------------------------------------------------------------------

			double *supply = new double[ORIG];
			double *demand = new double[DEST];
			double *c = new double[DEST*ORIG];	// linear costs origin-destination
			double *ub = new double[DEST*ORIG];	// upper bounds on the decision variables
			double total_supply;

			int k;
			double ran_x;
			double ran_y;

			//---------------------------------------------------------------------------------
			///////////////// DATA FILE BUILDING /////////////////////////////////////////////
			//---------------------------------------------------------------------------------

			if(SpaceDist == 1){

			//---------------------------------------------------------------------------------
			// Linear spatial distribution with linearly growing supplies and uniform demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
					supply[i] = ORIG*(1+i);
					total_supply = total_supply + supply[i];
					MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
					demand[j] = factor*(total_supply- eps)/DEST;
					MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-((int)round(ran_x*j)+1)) + abs(Strech*(i+1) - abs(((int)round(ran_y*j)+1)-(int)round(AA*DEST))))/Scal;
						ub[k] = demand[j]*c[k];

						MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
						k++;

					}
				}

			} else if(SpaceDist == 2){

			//---------------------------------------------------------------------------------
			// Linear spatial distribution with quadratically growing supplies and uniform demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
					supply[i] = ORIG*(1+i)*(1+i);
					total_supply = total_supply + supply[i];
					MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
					demand[j] = factor*(total_supply- eps)/DEST;
					MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-((int)round(ran_x*j)+1)) + abs(Strech*(i+1) - abs(((int)round(ran_y*j)+1)-(int)round(AA*DEST))))/Scal;
						ub[k] = demand[j]*c[k];

						MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
						k++;

					}
				}

			} else if(SpaceDist == 3){

			//---------------------------------------------------------------------------------
			// Circular distribution of customers with linearly growing supplies and uniform demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
					supply[i] = ORIG*(1+i);
					total_supply = total_supply + supply[i];
					MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
					demand[j] = factor*(total_supply- eps)/DEST;
					MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-(int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_x*j)/DEST)/2)) + abs(Strech*(i+1) - (int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_y*j)/DEST)/2)) )/(Scal);
						ub[k] = demand[j]*c[k];

						MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
						k++;

					}
				}

			} else if(SpaceDist == 4){

			//---------------------------------------------------------------------------------
			// Circular distribution of customers with quadratically growing supplies and uniform demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
					supply[i] = ORIG*(1+i)*(1+i);
					total_supply = total_supply + supply[i];
					MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
					demand[j] = factor*(total_supply- eps)/DEST;
					MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-(int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_x*j)/DEST)/2)) + abs(Strech*(i+1) - (int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_y*j)/DEST)/2)) )/(Scal);
						ub[k] = demand[j]*c[k];

						MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
						k++;

					}
				}

			} else if(SpaceDist == 5){

			//---------------------------------------------------------------------------------
			// Linear spatial distribution with linearly growing supplies and linearly growing demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
					supply[i] = ORIG*(1+i);
					total_supply = total_supply + supply[i];
					MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
					demand[DEST-j-1] = (factor*2*total_supply*(j + 1) - eps)/(DEST*(DEST+1));
					MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-((int)round(ran_x*j)+1)) + abs(Strech*(i+1) - abs(((int)round(ran_y*j)+1)-(int)round(AA*DEST))))/Scal;
						ub[k] = demand[j]*c[k];

						MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
						k++;

					}
				}

			} else if(SpaceDist == 6){

			//---------------------------------------------------------------------------------
			// Linear spatial distribution with quadratically growing supplies and linearly growing demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
					supply[i] = ORIG*(1+i)*(1+i);
					total_supply = total_supply + supply[i];
					MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
					demand[DEST-j-1] = (factor*2*total_supply*(j + 1)- eps)/(DEST*(DEST+1));
					MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-((int)round(ran_x*j)+1)) + abs(Strech*(i+1) - abs(((int)round(ran_y*j)+1)-(int)round(AA*DEST))))/Scal;
						ub[k] = demand[j]*c[k];

						MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
						k++;

					}
				}

			} else if(SpaceDist == 7){

			//---------------------------------------------------------------------------------
			// Circular distribution of customers with linearly growing supplies and linearly growing demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
					supply[i] = ORIG*(1+i);
					total_supply = total_supply + supply[i];
					MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
					demand[DEST-j-1] = (factor*2*total_supply*(j + 1) - eps)/(DEST*(DEST+1));
					MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-(int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_x*j)/DEST)/2)) + abs(Strech*(i+1) - (int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_y*j)/DEST)/2)) )/(Scal);
						ub[k] = demand[j]*c[k];

						MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
						k++;

					}
				}

			} else if(SpaceDist == 8){

			//---------------------------------------------------------------------------------
			// Circular distribution of customers with quadratically growing supplies and linearly growing demands
			//---------------------------------------------------------------------------------

				MyOutputFile << "c ---------------------------------- Node descriptor lines (supply)" << endl;

				total_supply = 0;
				for(int i=0; i<ORIG; i++){
					supply[i] = ORIG*(1+i)*(1+i);
					total_supply = total_supply + supply[i];
					MyOutputFile << "n1 " << i+1 << " " << supply[i] << endl;
				}

				MyOutputFile << "c ---------------------------------- Node descriptor lines (demand)" << endl;

				for(int j=0; j<DEST; j++){
					demand[DEST-j-1] = (factor*2*total_supply*(j + 1)- eps)/(DEST*(DEST+1));
					MyOutputFile << "n2 " << j+ORIG+1 << " " << -demand[j] << endl;
				}

				MyOutputFile << "c ---------------------------------- Arc descriptor lines (from, to, minflow, maxflow, cost)" << endl;

				k = 0;
				for(int i=0; i<ORIG; i++){
					ran_x = MM + (1-MM)*((double)My_Rand()/(4294967295));
					ran_y = MM + (1-MM)*((double)My_Rand()/(4294967295));

					for(int j=0; j<DEST; j++){

		        			c[k] = (abs(Strech*(i+1)-(int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_x*j)/DEST)/2)) + abs(Strech*(i+1) - (int)round(ORIG/2 + ORIG*cos((2*3.1415*ran_y*j)/DEST)/2)) )/(Scal);
						ub[k] = demand[j]*c[k];

						MyOutputFile << "a " << i+1 << " " << j+ORIG+1 << " " << 0 << " " << ub[k] << " " << c[k] << endl;
						k++;

					}
				}
			}

			delete[] c;
			delete[] ub;
			delete[] supply;
			delete[] demand; 


		}



	} // end try

	catch(...){
		// MyOutputFile << "Unknown exception" << endl;
	}

	MyOutputFile.close();

	return 0;

}

