#include "src/femwrks.h"

int main()
{
    cout << "-----------------------------------------" << endl;
	cout << "Running burger 3D case" << endl;
	cout << "-----------------------------------------" << endl;
    cout << "Constructing mesh..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(1,1,0.1);
    RESOLUTION d(20,20,2);
    MESH mesh(p, s, d, "H8");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials;
    materials.D = 0.1;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    bc.U.assignDBC(mesh, "x", 0, 0);
    bc.U.assignDBC(mesh, "x", 1, 0);
    bc.U.assignDBC(mesh, "y", 0, 0);
    bc.U.assignDBC(mesh, "y", 1, 0);

    cout << "Initializing field variables..."<< endl;
    POINT cntr(0.5,0.5,0);
    double a = 0.25; //radius of the implicit sphere
    vector<double> phi = setVector(mesh.TotalPoints, 1);
    for (int i = 1; i <= mesh.TotalPoints; i++)
    {phi[i] = minnum(phi[i],pow(mesh.points[i].x-cntr.x,2)+pow(mesh.points[i].y-cntr.y,2)+pow(mesh.points[i].z-cntr.z,2)-pow(a,2));}

    vector<double> U = setVector(mesh.TotalPoints, 0);
    for (int i = 1; i <= mesh.TotalPoints; i++)
    {
        if (phi[i] <= 0)
        {U[i] = 10;}
    }
    U = setVector(U, bc.U);

    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.SUPG = true;
    settings.stab_prmtr = 0.001;
    settings.transient = true;
    settings.dt = 0.0001;
    settings.nt = 10000;
    settings.prnt_freq = 10;
    settings.iters = 1;
    settings.CGiter = 1000;
    settings.tol = 0.0001;

    cout << "Running solver..."<< endl;
    //solveBurger(mesh, materials, bc, settings, U);
    solveBurgerEBE(mesh, materials, bc, settings, U);

    return 0;
}
