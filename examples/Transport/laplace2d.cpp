#include "src/femwrks.h"

int main()
{
    cout << "-----------------------------------------" << endl;
	cout << "2D laplace equation" << endl;
	cout << "-----------------------------------------" << endl;
	
    cout << "Constructing mesh..."<< endl;
    //MESH mesh("msh/PlateWithHole2D/mesh_400.dat");
    MESH mesh("msh/PlateWithHole2D_O2/mesh_1200.dat");
	
	cout << "Assigning material properties..."<< endl;
    MATERIALS materials;
    materials.D = 1;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITION bc;
    //bc.assignDBC(mesh, "msh/PlateWithHole2D/top_400.dat",0);
    //bc.assignDBC(mesh, "msh/PlateWithHole2D/bottom_400.dat",10);
    
	bc.assignDBC(mesh, "msh/PlateWithHole2D_O2/top_1200.dat",0);
    bc.assignDBC(mesh, "msh/PlateWithHole2D_O2/bottom_1200.dat",10);
    
	cout << "Initializing field variables..."<< endl;
    vector<double> Vx = setVector(mesh.TotalPoints, 0);
    vector<double> Vy = setVector(mesh.TotalPoints, 0);
    vector<double> Vz = setVector(mesh.TotalPoints, 0);
    vector<double> F = setVector(mesh.TotalPoints, 0);
    vector<double> U = setVector(mesh.TotalPoints, 0);

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
    clock_t start, end;
    start = clock();
    //solveTransport(mesh, materials, bc, settings, U, Vx, Vy, Vz, F);
    solveTransportEBE(mesh, materials, bc, settings, U, Vx, Vy, Vz, F);
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by program is : " << time_taken << setprecision(10) << " seconds" << endl;
}
