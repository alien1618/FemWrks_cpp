#ifndef FEMWRKS_H
#define FEMWRKS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <omp.h>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <tuple>

using namespace std;

struct SOLVER_SETTINGS
{
    double tol;			//convergence tolerance
    int CGiter;  		//conjugate gradient itererations
    int iters;   		//main solver iterations
    bool SUPG; 			//streamline upwind petrov galerkin
    double stab_prmtr; 	//supg parameter
    bool transient; 	//transient or steady-state
    double dt;  		//time step
    int nt;  			//total time steps
    int prnt_freq;  	//print solution frequency
};

struct MATERIALS
{
    double D; 			//diffusivity = (K/(ro*Cp))
    double ro; 			//density
    double Nu; 			//viscosity
    double YTS; 		//Yeild tensile strength
    double MOE; 		//modulous of elasticity
    double nu; 			//poisson's ratio
    double thickness; 	//material thickness
};

class POINT
{
private:

public:
    double x;
    double y;
    double z;
    
    POINT();
    POINT(double);
    POINT(double, double);
    POINT(double, double, double);
    double distance(POINT point);
    ~POINT();
};
class DOMAIN
{
private:

public:
    double x;
    double y;
    double z;
    
	DOMAIN();
	DOMAIN(double);
    DOMAIN(double, double);
    DOMAIN(double, double, double);
    ~DOMAIN();
};
class RESOLUTION
{
private:

public:
    int x;
    int y;
    int z;
    
    RESOLUTION();
    RESOLUTION(int);
    RESOLUTION(int, int);
    RESOLUTION(int, int, int);
    ~RESOLUTION();
};

class MESH
{

private:

public:
    MESH ();
    MESH (string input_txt);
    MESH (POINT p, DOMAIN l, RESOLUTION s, string elemtype);

    //----------------------------------------------------------------
    // MESH DATA
    //----------------------------------------------------------------
    int dim;
    int ElementOrder;
    
    vector<POINT> points;
    int TotalPoints;

    vector<vector<int> > lines;    
    int LineNds;
    int TotalLines;
    
    vector<vector<int> > surfaces;
    int SurfaceShape;
    int SurfaceNds;
    int TotalSurfaces;
    
    vector<vector<int> > volumes;
    int VolumeShape;
    int VolumeNds;
    int TotalVolumes;

    //------------------------------------------
    //Generation of Structured Grids Functions
    //------------------------------------------
    void generateGridNodes1D(POINT p, DOMAIN l, RESOLUTION div);
    void generateGridNodes2D(POINT p, DOMAIN l, RESOLUTION div);
    void generateGridNodes3D(POINT p, DOMAIN l, RESOLUTION div);
    void generateGridQuad(RESOLUTION div);
    void generateGridTri(RESOLUTION div);
    void generateGridTri2(RESOLUTION div);
    void generateGridTet(POINT p, DOMAIN l, RESOLUTION div);
    void generateGridHex(POINT p, DOMAIN l, RESOLUTION div);

    //------------------------------------------
    //Printing Mesh Data
    //------------------------------------------
    void printMeshVTK(vector<double> U, string filename, int t);
    void printMeshVTK(vector<int> U, string filename, int t);
    void printVectorsVTK(vector<double> Vx, vector<double> Vy, vector<double> Vz, string filename, int t);
    void printTXT(vector<double> U, string filename, int t);
	void print();
    ~MESH ();
};

class GAUSS
{
private:

public:
    int element_integpoints;
    vector<double> eps;
    vector<double> eta;
    vector<double> zta;
    vector<double> w;
    void tetra(int num_pts);
    void hexa(int num_pts);
    void quad(int num_pts);
    void tri(int num_pts);
    void line(int num_pts);
    void generatePnts(int dim, int elem_type, int level);

    ~GAUSS ();
};

class KERNEL
{
private:
    void jacobianVolume(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, vector<double> dNdeps, vector<double> dNdeta, vector<double> dNdzta);
    void jacobianSurface(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, vector<double> dNdeps, vector<double> dNdeta);

public:
    int elem_num;  
    vector<double> N;
    vector<double> dNdx;
    vector<double> dNdy;
    vector<double> dNdz;
	vector<vector<double> > dNdNx;
	vector<vector<double> > dNdNy;
	vector<vector<double> > dNdNz;
	vector<vector<double> > NdNx;
	vector<vector<double> > NdNy;
	vector<vector<double> > NdNz;
	vector<vector<double> > NN;
	vector<vector<double> > MK;

	double quad_weight;
    double det_jac;

    POINT pnt;
    void computeLagrange(int dim, vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps, double eta, double zta);
    void T4(vector<POINT> points, vector<int> nbrnds, int TotNbrNds,  double eps, double eta, double zta);
    void T10(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps, double eta, double zta);
    void H8(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps, double eta, double zta);
    void H20(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps, double eta, double zta);
    void Q4(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps,double eta);
    void Q8(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps,double eta);
    void T3(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps,double eta);
    void T6(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps,double eta);
    void lineO1(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps);

    ~KERNEL ();
};
vector<KERNEL> computeInterpolants(int dim, vector<POINT> points, vector<vector<int> > elements, int totelems, int elemnds, GAUSS gauss);

class BOUNDARY_CONDITION
{
private:

public:
    BOUNDARY_CONDITION();
    
    int total;
	vector<int> points;
	vector<double> values;
    void assignDBC(MESH mesh, string salomeface_fname, double ValuePlusPhi);
    void assignDBC(MESH mesh, string edge, double location, double value);

    ~BOUNDARY_CONDITION ();
};

class BOUNDARY_CONDITIONS
{
private:

public:
    BOUNDARY_CONDITION Vx;
    BOUNDARY_CONDITION Vy;
    BOUNDARY_CONDITION Vz;
    BOUNDARY_CONDITION P;
    BOUNDARY_CONDITION U;
    BOUNDARY_CONDITION Fx;
    BOUNDARY_CONDITION Fy;
    BOUNDARY_CONDITION Fz;
};

vector<double> setVector(vector<double> U, BOUNDARY_CONDITION BC);
vector<double> setVector(int TotalPoints, double val);

class GLOBAL_MATRICES
{
private:

public:

    vector<vector<double> > K;
    vector<vector<double> > M;
    vector<double> F;

    void assembleDiffusionMatrix(bool SUPG, int TotalMeshPoints, vector<vector<int>> elements, int elem_nds, vector<KERNEL> interps, int gauss_TotalPoints, double D, vector<double> Vx1, vector<double> Vy1, vector<double> Vz1, double stab_prmtr);
    void assembleMassMatrix(bool SUPG, int TotalMeshPoints, vector<vector<int>> elements, int elem_nds, vector<KERNEL> interps, int gauss_TotalPoints, vector<double> Vx1, vector<double> Vy1, vector<double> Vz1, double stab_prmtr);
    void assembleStiffnessMatrix(int dim, int load_type, double E, double nu, double t, int TotalNds, vector<KERNEL> kernels, int totquadpnts,  vector<vector<int> > elements, int elemnds);
    void assembleForceMatrix(bool SUPG, int TotalMeshPoints, vector<vector<int>> elements, int elem_nds, vector<KERNEL> interps, int gauss_TotalPoints, vector<double> Vx1, vector<double> Vy1, vector<double> Vz1, vector<double> Q1, double stab_prmtr);

    void applyDBCDirectSteadyState(int DOF, int TotalPoints, vector<int> unique_points, vector<double> unique_values, int TotalUniquePoints);
    void applyDBCDirectTransient(int TotalPoints, vector<int> unique_points, vector<double> unique_values, int TotalUniquePoints);
    void applyDBCPenaltyTransient(int TotalPoints, vector<int> unique_points, vector<double> unique_values, int TotalUniquePoints);

    ~GLOBAL_MATRICES();
};
vector<double> advanceInTime(int transient, vector<double> U_now, GLOBAL_MATRICES global_matrices, int TotalPoints, double dt, int CGiter,  double tol);


//----------------------------------------------------------------------
// MISC FUNCTIONS
//----------------------------------------------------------------------
void printElements(vector<vector<int> > surfaces, int SurfaceNds, int TotalSurfaces, string filename);
void pltctrl(vector<POINT> points, int TotalPoints, int n_t, int print_frequency);

double mod(double number, double divnum);
double maxnum(double a, double b);
double minnum(double a, double b);

void randomSeed();
double random0to1();

template <class type>
type maximum(vector<type> V, int z)
{
    type max = V[1];
        #pragma omp parallel for
        for (int i = 2; i <= z; i++)
        {
            type number1 = V[i];
            if (number1 >= max)
            {max = number1;}
        }
    return max;
}

template <typename type>
type minimum(vector<type> V, int z)
{
    type min = V[1];
        #pragma omp parallel for
        for (int i = 2; i <= z; i++)
        {
            type number1 = V[i];
            if (number1 <= min)
            {min = number1;}
        }
    return min;
}

double det_2X2(vector<vector<double> > mat);
double det_3X3(vector<vector<double> > mat);
vector<vector<double> > inverse_2X2(vector<vector<double> > A);
vector<vector<double> > inverse_3X3(vector<vector<double> > A);

vector<double> add(vector<double> A, vector<double> B, int size);
vector<vector<double> > add(vector<vector<double> > A, vector<vector<double> > B, int clmns, int rows);

vector<vector<double> > multiply(vector<vector<double> > A, int A_clmns, int A_rows, vector<vector<double> > B, int B_clmns, int B_rows);
vector<double> multiply(vector<vector<double> > A, int A_Xclmns, int A_rows, vector<double> B, int B_size);
vector<vector<double> > multiply_transp(vector<double> A, vector<double> B, int size);
vector<double> multiply(vector<double> A, int A_size, double number);
vector<vector<double> > multiply(vector<vector<double> > A, int A_Xclmns, int A_rows,  double number);
vector<vector<double> > transpose(vector<vector<double> > A, int A_Xclmns, int A_Yrows);

vector<double> solveTransientImplicitPCG(vector<vector<double> > K, vector<vector<double> > N, vector<double> F, vector<double> conc_now, int TotalPoints, double dt, int CGiter, double tol);
vector<double> solveTransientExplicitPCG(vector<vector<double> > K, vector<vector<double> > N, vector<double> F, vector<double> conc_now, int TotalPoints, double dt, int CGiter, double tol);
vector<double> solvePCG(vector<vector<double> > A, vector<double> b, int mat_size, int CGiter, double tol);
vector<vector<double> > inv(vector<vector<double> > a, int size);
vector<double> solveInverse(vector<vector<double> > K, vector<double> F, int TotalPoints);

//----------------------------------------------------------------------
// GLOBAL-MATRIX FREE OPERATIONS
//----------------------------------------------------------------------
struct NBR
{
    vector<int> nbrs;
    int TotalNbrs;
};
vector<NBR> computeNodeKernels(int TotalPoints, vector<vector<int> > surfaces, int TotalSurfaces, int SurfaceNds);

tuple <vector<int>, int, vector<int>, int> collectGaussPointsWithBC(vector<int> bc_ndnum, int BC_TotNds, vector<vector<int>> elements, int SurfaceNds, vector<KERNEL> kernels, int TotalGaussPoints);

vector<double> solveElementByElement(vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, double D, int TotNds, vector<KERNEL> kernels, int TotalGaussPoints, vector<vector<int>> elements, int SurfaceNds, vector<NBR> ndkernels, double dt, int transient, bool SUPG, double stab_prmtr, int CGiter, double conv_tol, vector<int> bc_ndnum, vector<double> bc_ndval, int BC_TotNds, vector<int> g_bc, int tot_g_bc, vector<int> g_bc2, int tot_g_bc2);

vector<double> solveElementByElement(int dim, int DOF, double E, double nu, double t, int TotNds, vector<KERNEL> kernels, int TotalGaussPoints, vector<vector<int>> elements, int SurfaceNds, vector<NBR> ndkernels, int CGiter, double conv_tol, vector<double> F, vector<int> bc_ndnum, vector<double> bc_ndval, int BC_TotNds, vector<int> g_bc, int tot_g_bc, vector<int> g_bc2, int tot_g_bc2);

vector<double> solvePCG(int DOF, vector<double> F, vector<NBR> ndkernels, int TotNds, vector<vector<int>> elements, int SurfaceNds, vector<KERNEL> kernels, vector<int> bc_ndnum, vector<double> bc_ndval, int BC_TotNds, vector<int> g_bc, int tot_g_bc, vector<int> g_bc2, int tot_g_bc2, int CGiter, double conv_tol);


//----------------------------------------------------------------------
// SOLVER FUNCTIONS
//----------------------------------------------------------------------
vector<vector<double> > computeD2D(int load_type, double E, double nu, double t);
vector<vector<double> > computeD3D(double E, double nu);
vector<vector<double> > computeB2D(int SurfaceNds, vector<double> dNdx, vector<double> dNdy);
vector<vector<double> > computeB3D(int VolumeNds, vector<double> dNdx, vector<double> dNdy, vector<double> dNdz);

void solveElasticity(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings);
void solveElasticityEBE(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings);

void solveTransport(MESH mesh, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q);
void solveTransportEBE(MESH mesh, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q);

void solveBurger(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> Vx);
void solveBurgerEBE(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> Vx);

void solveNavierStokes(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC,  SOLVER_SETTINGS settings, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Qx, vector<double> Qy, vector<double> Qz);
void solveNavierStokesEBE(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Qx, vector<double> Qy, vector<double> Qz);

tuple<vector<double>, vector<double>, vector<double>, vector<double> > computePressure(BOUNDARY_CONDITION BC_P, double dt, vector<double> Vx_inter, vector<double> Vy_inter, vector<double> Vz_inter, int TotalPoints, vector<vector<int>> elements, int elem_nds, vector<KERNEL> interpolants, int TotalGaussPoints, double density, int CGiter, double conv_tol);
tuple<vector<double>, vector<double>, vector<double>, vector<double> > computePressureEBE(vector<double> Vx_inter, vector<double> Vy_inter, vector<double> Vz_inter, int TotNds, vector<KERNEL> kernels, int gauss_TotalPoints, vector<vector<int>> elements, int ElementNds, vector<NBR> ndkernels, double ro, double dt, int CGiter, double conv_tol, vector<int> bc_ndnum, vector<double> bc_ndval, int BC_TotNds, vector<int> g_bc, int tot_g_bc, vector<int> g_bc2, int tot_g_bc2);

#endif
