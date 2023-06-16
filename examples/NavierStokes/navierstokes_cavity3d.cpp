#include "src/femwrks.h"

int main()
{
    cout << "-------------------------------------" << endl;
	cout << "3D cavity flow" << endl;
	cout << "-------------------------------------" << endl;
	cout << "Constructing mesh..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(1,1,0.1);
    RESOLUTION d(20,20,2);
    MESH mesh(p, s, d, "T4");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials;
    materials.Nu = 0.1;
    materials.ro = 1;

    cout << "Assigning boundary conditions..."<< endl;
    double vel = 10;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(mesh, "x", 0, 0);
    bc.Vy.assignDBC(mesh, "x", 0, 0);
    bc.Vz.assignDBC(mesh, "x", 0, 0);
    bc.Vx.assignDBC(mesh, "x", 1, 0);
    bc.Vy.assignDBC(mesh, "x", 1, 0);
    bc.Vz.assignDBC(mesh, "x", 1, 0);
    bc.Vx.assignDBC(mesh, "y", 0, 0);
    bc.Vy.assignDBC(mesh, "y", 0, 0);
    bc.Vz.assignDBC(mesh, "y", 0, 0);
    bc.Vx.assignDBC(mesh, "y", 1, vel);
    bc.Vy.assignDBC(mesh, "y", 1, 0);
    bc.Vz.assignDBC(mesh, "y", 1, 0);
    bc.P.assignDBC(mesh, "y", 0, 0);

    cout << "Initializing field variables..."<< endl;
    vector<double> Vx = setVector(mesh.TotalPoints, 1e-6);
    vector<double> Vy = setVector(mesh.TotalPoints, 1e-6);
    vector<double> Vz = setVector(mesh.TotalPoints, 1e-6);
    vector<double> Fx = setVector(mesh.TotalPoints, 0);
    vector<double> Fy = setVector(mesh.TotalPoints, 0); 
    vector<double> Fz = setVector(mesh.TotalPoints, 0);  
    
    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.SUPG = true;
    settings.stab_prmtr = 0.01;
    settings.CGiter = 1000;
    settings.tol = 0.001;
    settings.transient = true;
    settings.dt = 0.0005;
    settings.nt = 10000;
    settings.prnt_freq = 10;

    cout << "Running solver..."<< endl;
    //solveNavierStokes(mesh, materials, bc, settings, Vx, Vy, Vz, Fx, Fy, Fz);
    solveNavierStokesEBE(mesh, materials, bc, settings, Vx, Vy, Vz, Fx, Fy, Fz);
    
    return 0;
}
