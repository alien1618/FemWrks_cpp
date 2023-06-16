#include "src/femwrks.h"

int main()
{
    cout << "-----------------------------------------" << endl;
	cout << "2D flow past cylinder" << endl;
	cout << "-----------------------------------------" << endl;
	cout << "Constructing mesh..."<< endl;
    MESH mesh("mesh/FlowPastCylinder2D/mesh_627.dat");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials;
    materials.Nu = 0.1;
    materials.ro = 1;

    cout << "Assigning boundary conditions..."<< endl;
    double Vxb = 1;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(mesh,"mesh/FlowPastCylinder2D/inlet_627.dat", Vxb);
    bc.Vy.assignDBC(mesh,"mesh/FlowPastCylinder2D/inlet_627.dat", 0);
    bc.Vz.assignDBC(mesh,"mesh/FlowPastCylinder2D/inlet_627.dat", 0);
    bc.Vx.assignDBC(mesh,"mesh/FlowPastCylinder2D/walls_627.dat", 0);
    bc.Vy.assignDBC(mesh,"mesh/FlowPastCylinder2D/walls_627.dat", 0);
    bc.Vz.assignDBC(mesh,"mesh/FlowPastCylinder2D/walls_627.dat", 0);
    bc.P.assignDBC(mesh,"mesh/FlowPastCylinder2D/outlet_627.dat", 0);

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
    settings.stab_prmtr = 0.001;
    settings.CGiter = 1000;
    settings.tol = 0.001;
    settings.iters = 1;
    settings.transient = true;
    settings.dt = 0.0001;
    settings.nt = 50000;
    settings.prnt_freq = 1;
    
    cout << "Running solver..."<< endl;
    solveNavierStokes(mesh, materials, bc, settings, Vx, Vy, Vz, Fx, Fy, Fz);

    return 0;
}
