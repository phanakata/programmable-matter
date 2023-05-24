#define PERIOD 10000
#define NMAX 200000
#define MAXPARTICLETYPE 10
#define a 1.0

extern int NX, NY, N,Nb,Nd,bondGroup[NMAX*2],dihedralGroup[NMAX*4],CLONE, CALC;
extern float position[NMAX*3];
extern uint32_t particleID[NMAX];
extern char particleType[MAXPARTICLETYPE][2];
extern double fmesh[NMAX], DELTA;
extern int countfmesh[NMAX];

//extern double bendingEner[NMAX];
//extern double bondHarmonicEner[NMAX];
//extern double total_DHE,total_BHE;
//extern double hgt_fluctuation[NMAX];
//extern double h_avg_node[NMAX];
//extern double hgt_fluctuation[NMAX];
//extern double h_width[FRAMES/2][2*NX];
