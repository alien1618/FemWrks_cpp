#include "src/femwrks.h"

int main()
{
    cout << "-----------------------------------------" << endl;
	cout << "2D advection-diffusion" << endl;
	cout << "-----------------------------------------" << endl;
	cout << "Constructing mesh..."<< endl;
    POINT p(-1,-1);
    DOMAIN s(2,2);
    RESOLUTION d(50,50);
    MESH mesh(p, s, d, "T3");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials;
    materials.D = 0.0;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITION bc;
    //bc.assignDBC(mesh, "x", 0, 1);
    //bc.assignDBC(mesh, "y", 0, 1);
    //bc.assignDBC(mesh, "x", 1, 1);
    //bc.assignDBC(mesh, "y", 1, 1);

    cout << "Initializing field variables..."<< endl;
    POINT cntr(0,-0.5);
    double a = 0.25;  
    vector<double> phi = setVector(mesh.TotalPoints, 1);
    for (int i = 1; i <= mesh.TotalPoints; i++)
    {phi[i] = minnum(phi[i],pow(mesh.points[i].x-cntr.x,2)+pow(mesh.points[i].y-cntr.y,2)+pow(mesh.points[i].z-cntr.z,2)-pow(a,2));}
    
    vector<double> Vx = setVector(mesh.TotalPoints, 0);
    vector<double> Vy = setVector(mesh.TotalPoints, 0);
    vector<double> Vz = setVector(mesh.TotalPoints, 0);
    vector<double> F = setVector(mesh.TotalPoints, 0);
    vector<double> U = setVector(mesh.TotalPoints, 0);
    double u = 0.3;
    double v = 0.3;
    for (int i = 1; i <= mesh.TotalPoints; i++)
    {
        //simple rotation around center
        Vx[i] = u*(mesh.points[i].y);
        Vy[i] = -v*(mesh.points[i].x);
        if (phi[i] <= 0)
        {U[i] = 1;}
    }
    mesh.printMeshVTK(U,"U", 0);

    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.transient = true;
    settings.dt = 0.01;
    settings.nt = 1000;
    settings.SUPG = true;
    settings.stab_prmtr = 0.05;
    settings.prnt_freq = 10;
    settings.iters = 5;
    settings.CGiter = 1000;
    settings.tol = 0.0001;

    cout << "Running solver..."<< endl;
    //solveTransport(mesh, materials, bc, settings, U, Vx, Vy, Vz, F);
    solveTransportEBE(mesh, materials, bc, settings, U, Vx, Vy, Vz, F);

    return 0;
}
