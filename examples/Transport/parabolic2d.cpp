#include "src/femwrks.h"

int main()
{
    cout << "-----------------------------------------" << endl;
	cout << "2D parabolic equation" << endl;
	cout << "-----------------------------------------" << endl;
    
	cout << "Constructing mesh..."<< endl;
    MESH mesh("msh/PlateWithHole2D/mesh_400.dat");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials;
    materials.D = 0.1;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITION bc;
    bc.assignDBC(mesh, "msh/PlateWithHole2D/top_400.dat",0);
    bc.assignDBC(mesh, "msh/PlateWithHole2D/bottom_400.dat",10);

    cout << "Initialize distribution of field variables in the domain..."<< endl;
    vector<double> Vx = setVector(mesh.TotalPoints, 0);
    vector<double> Vy = setVector(mesh.TotalPoints, 0);
    vector<double> Vz = setVector(mesh.TotalPoints, 0);
    vector<double> F = setVector(mesh.TotalPoints, 0);
    vector<double> U = setVector(mesh.TotalPoints, 0);
        
    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.transient = true;
    settings.dt = 0.01;
    settings.nt = 100;
    settings.prnt_freq = 1;
    settings.SUPG = false;
    settings.stab_prmtr = 0.0;
    settings.iters = 1;
    settings.CGiter = 1000;
    settings.tol = 0.0001;

    cout << "Running solver..."<< endl;
    //solveTransport(mesh, materials, bc, settings, U, Vx, Vy, Vz, F);
    solveTransportEBE(mesh, materials, bc, settings, U, Vx, Vy, Vz, F);

    return 0;
}
