#include "src/femwrks.h"

int main()
{
    cout << "-------------------------------------" << endl;
	cout << "2D deformation-stress analysis" << endl;
	cout << "-------------------------------------" << endl;
	cout << "Constructing mesh..."<< endl;
    //MESH mesh("msh/PlateWithHole2D/mesh_400.dat");
    MESH mesh("msh/PlateWithHole2D_O2/mesh_1200.dat");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials;
    materials.YTS = 150;
    materials.MOE = 10000;
    materials.nu = 0.3;
    materials.thickness = 1;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    /*
    bc.U.assignDBC(mesh, "msh/PlateWithHole2D/bottom_400.dat",0);
    bc.Fx.assignDBC(mesh, "msh/PlateWithHole2D/top_400.dat", 0);
    bc.Fy.assignDBC(mesh, "msh/PlateWithHole2D/top_400.dat", -50);
    */
    
    bc.U.assignDBC(mesh, "msh/PlateWithHole2D_O2/bottom_1200.dat",0);
    bc.Fx.assignDBC(mesh, "msh/PlateWithHole2D_O2/top_1200.dat", 0);
    bc.Fy.assignDBC(mesh, "msh/PlateWithHole2D_O2/top_1200.dat", -20);
    
    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.transient = false;
    settings.dt = 1;
    settings.nt = 1;
    settings.iters = 1;
    settings.CGiter = 1000;
    settings.tol = 0.0001;
    settings.SUPG = false;
    settings.stab_prmtr = 0;
    settings.prnt_freq = 1;

    cout << "Running solver..."<< endl;
    clock_t start, end;
    start = clock();
    solveElasticity(mesh, materials, bc, settings);
    //solveElasticityEBE(mesh, materials, bc, settings);
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by program is : " << time_taken << setprecision(10) << " seconds" << endl;
    
    return 0;
}
