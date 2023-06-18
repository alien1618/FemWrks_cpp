#include "src/femwrks.h"

int main()
{
    cout << "-----------------------------------------" << endl;
	cout << "3D poisson equation" << endl;
	cout << "-----------------------------------------" << endl;
	cout << "Constructing mesh..."<< endl;
    POINT p(-1,-1,0);
    DOMAIN s(2,2,0.1);
    RESOLUTION d(20, 20, 2);
    MESH mesh(p, s, d, "T4");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials;
    materials.D = 1;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITION bc;
    bc.assignDBC(mesh, "x", -1, 1);
    bc.assignDBC(mesh, "y", -1, 1);
    bc.assignDBC(mesh, "x", 1, 1);
    bc.assignDBC(mesh, "y", 1, 1);

    cout << "Initializing field variables..."<< endl;
    vector<double> Vx = setVector(mesh.TotalPoints, 0);
    vector<double> Vy = setVector(mesh.TotalPoints, 0);
    vector<double> Vz = setVector(mesh.TotalPoints, 0);
    vector<double> F = setVector(mesh.TotalPoints, 0);
    vector<double> U = setVector(mesh.TotalPoints, 0);
    for (int i = 1; i <= mesh.TotalPoints/2; i++)
    {F[i] = -10;} 
 
    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.transient = false;
    settings.dt = 1;
    settings.nt = 1;
    settings.prnt_freq = 1;
    settings.SUPG = false;
    settings.stab_prmtr = 0.0;
    settings.iters = 1;
    settings.CGiter = 1000;
    settings.tol = 0.0001;

    cout << "Running solver..."<< endl;
    solveTransportEBE(mesh, materials, bc, settings, U, Vx, Vy, Vz, F);

    return 0;
}
