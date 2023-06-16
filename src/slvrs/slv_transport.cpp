#include "../femwrks.h"

void solveTransport(MESH mesh, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q)
{
    cout << "SOLVER START..." << endl;
    pltctrl(mesh.points, mesh.TotalPoints, settings.nt, settings.prnt_freq);
    if (mesh.dim == 2)
    {printElements(mesh.surfaces, mesh.SurfaceNds, mesh.TotalSurfaces, "elements");}
    
    cout << "Computing quadrature points..." << endl;
    GAUSS gauss;
    int TotQuadPnts = 0;
    switch(mesh.dim)
    {
        case 2:
        {
            gauss.generatePnts(mesh.dim, mesh.SurfaceShape, mesh.ElementOrder);
            TotQuadPnts = gauss.element_integpoints*mesh.TotalSurfaces;
        }break;
        case 3:
        {
            gauss.generatePnts(mesh.dim, mesh.VolumeShape, mesh.ElementOrder);
            TotQuadPnts = gauss.element_integpoints*mesh.TotalVolumes;
        }break;
    }
    
    cout << "Computing interpolants..." << endl;
    vector<KERNEL> kernels;
    switch(mesh.dim)
    {
        case 2:
        {kernels = computeInterpolants(mesh.dim, mesh.points, mesh.surfaces, mesh.TotalSurfaces, mesh.SurfaceNds, gauss);}break;
        case 3:
        {kernels = computeInterpolants(mesh.dim, mesh.points, mesh.volumes, mesh.TotalVolumes, mesh.VolumeNds, gauss);}break;
    }

    cout << "Computing global matrices..." << endl;
    GLOBAL_MATRICES global_matrices;
    switch(mesh.dim)
    {
        case 2:
        {
           global_matrices.assembleDiffusionMatrix(settings.SUPG, mesh.TotalPoints, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts, materials.D, Vx, Vy, Vz, settings.stab_prmtr);
           global_matrices.assembleForceMatrix(settings.SUPG, mesh.TotalPoints, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts, Vx, Vy, Vz, Q, settings.stab_prmtr);
        }break;
        case 3:
        {
            global_matrices.assembleDiffusionMatrix(settings.SUPG, mesh.TotalPoints, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts, materials.D, Vx, Vy, Vz, settings.stab_prmtr);
            global_matrices.assembleForceMatrix(settings.SUPG, mesh.TotalPoints, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts, Vx, Vy, Vz, Q, settings.stab_prmtr);
        }break;
    }
    
    cout << "Applying boundary conditions to global matrices..." << endl;
    if (settings.transient == true)
    {
        switch(mesh.dim)
        {
            case 2:
            {
                global_matrices.assembleMassMatrix(settings.SUPG, mesh.TotalPoints, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts, Vx, Vy, Vz, settings.stab_prmtr);
            }break;
            case 3:
            {
                global_matrices.assembleMassMatrix(settings.SUPG, mesh.TotalPoints, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts, Vx, Vy, Vz, settings.stab_prmtr);
            }break;
            global_matrices.applyDBCDirectTransient(mesh.TotalPoints, BC.points, BC.values, BC.total);
        }
    }
    else
    {global_matrices.applyDBCDirectSteadyState(1, mesh.TotalPoints, BC.points, BC.values, BC.total);}

    U = setVector(U, BC);
    mesh.printMeshVTK(U,"U", 0);
    if (settings.transient == false)
    {
        cout << "Computing steady-state solution..." << endl;
        U = advanceInTime(settings.transient, U, global_matrices, mesh.TotalPoints, settings.dt, settings.CGiter, settings.tol);
        mesh.printMeshVTK(U, "U", 1);
        mesh.printTXT(U, "U", 1);
    }
    else
    {
        int cntr = 1;
        for (int t = 1; t <= settings.nt; t++)
        {
            cout << "Processing time step " << t << " of " << settings.nt <<  endl;

            cout << "Solving the system of equations..." << endl;
            U = advanceInTime(1, U, global_matrices, mesh.TotalPoints, settings.dt, settings.CGiter, 0.0001);
            if (t/settings.prnt_freq == cntr)
            {
                cntr++;
                mesh.printMeshVTK(U, "U", t);
                mesh.printTXT(U, "U", t);
            }
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}
void solveTransportEBE(MESH mesh, MATERIALS materials, BOUNDARY_CONDITION BC, SOLVER_SETTINGS settings, vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q)
{
    cout << "SOLVER START..." << endl;
    pltctrl(mesh.points, mesh.TotalPoints, settings.nt, settings.prnt_freq);
    if (mesh.dim == 2)
    {printElements(mesh.surfaces, mesh.SurfaceNds, mesh.TotalSurfaces, "elements");}
    
	cout << "Collect element support for each node..." << endl;
    vector<NBR> ndkernels;
    switch(mesh.dim)
    {
        case 2:
	    {ndkernels = computeNodeKernels(mesh.TotalPoints, mesh.surfaces, mesh.TotalSurfaces, mesh.SurfaceNds);}break;
        case 3:
	    {ndkernels = computeNodeKernels(mesh.TotalPoints, mesh.volumes, mesh.TotalVolumes, mesh.VolumeNds);}break;
    }
	cout << "Computing quadrature points..." << endl;
    GAUSS gauss;
    int TotQuadPnts = 0;
    switch (mesh.dim)
    {
        case 2:
        {
            gauss.generatePnts(mesh.dim, mesh.SurfaceShape, mesh.ElementOrder);
            TotQuadPnts = gauss.element_integpoints*mesh.TotalSurfaces;
        }break;
        case 3:
        {
            gauss.generatePnts(mesh.dim, mesh.VolumeShape, mesh.ElementOrder);
            TotQuadPnts = gauss.element_integpoints*mesh.TotalVolumes;
        }break;
    }
    
    cout << "Computing interpolants..." << endl;
    vector<KERNEL> kernels;
    switch(mesh.dim)
    {
        case 2:
        {kernels = computeInterpolants(mesh.dim, mesh.points, mesh.surfaces, mesh.TotalSurfaces, mesh.SurfaceNds, gauss);}break;
        case 3:
        {kernels = computeInterpolants(mesh.dim, mesh.points, mesh.volumes, mesh.TotalVolumes, mesh.VolumeNds, gauss);}break;
    }
    
	cout << "Collecting all boundary nodes and gauss points with boundary conditions..." << endl;
    vector<int>g_bc; int tot_g_bc=0; vector<int> g_bc2; int tot_g_bc2=0;
    switch(mesh.dim)
    {
        case 2:
        {tie(g_bc, tot_g_bc, g_bc2, tot_g_bc2) = collectGaussPointsWithBC(BC.points, BC.total, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts);}break;
        case 3:
        {tie(g_bc, tot_g_bc, g_bc2, tot_g_bc2) = collectGaussPointsWithBC(BC.points, BC.total, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts);}break;
    }

    cout << "Looping over time..." << endl;
	int cntr = 1;
	mesh.printMeshVTK(U, "U", 0);
    mesh.printTXT(U, "U", 0);
    for (int t = 1; t <= settings.nt; t++)
    {
        cout << "Processing time step " << t << " of " << settings.nt <<  endl;
        switch(mesh.dim)
        {
            case 2:
            {U = solveElementByElement(U, Vx, Vy, Vz, Q, materials.D, mesh.TotalPoints, kernels, TotQuadPnts, mesh.surfaces, mesh.SurfaceNds, ndkernels, settings.dt, settings.transient, settings.SUPG, settings.stab_prmtr, settings.CGiter, settings.tol, BC.points, BC.values, BC.total, g_bc, tot_g_bc, g_bc2, tot_g_bc2);}break;
            case 3:
            {U = solveElementByElement(U, Vx, Vy, Vz, Q, materials.D, mesh.TotalPoints, kernels, TotQuadPnts, mesh.volumes, mesh.VolumeNds, ndkernels, settings.dt, settings.transient, settings.SUPG, settings.stab_prmtr, settings.CGiter, settings.tol, BC.points, BC.values, BC.total, g_bc, tot_g_bc, g_bc2, tot_g_bc2);}break;
        }
        U = setVector(U, BC);
        if (t/settings.prnt_freq == cntr)
        {
            cntr++;
            mesh.printMeshVTK(U, "U", t);
            mesh.printTXT(U, "U", t);
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}

void solveBurger(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> U)
{
    cout << "SOLVER START..." << endl;
    pltctrl(mesh.points, mesh.TotalPoints, settings.nt, settings.prnt_freq);
    if (mesh.dim == 2)
    {printElements(mesh.surfaces, mesh.SurfaceNds, mesh.TotalSurfaces, "elements");}
    vector<double> zero(mesh.TotalPoints,0);   
    
	cout << "Computing quadrature points..." << endl;
    GAUSS gauss;
    int TotQuadPnts = 0;
    switch(mesh.dim)
    {
        case 2:
        {
            gauss.generatePnts(mesh.dim, mesh.SurfaceShape, mesh.ElementOrder);
            TotQuadPnts = gauss.element_integpoints*mesh.TotalSurfaces;
        }break;
        case 3:
        {
            gauss.generatePnts(mesh.dim, mesh.VolumeShape, mesh.ElementOrder);
            TotQuadPnts = gauss.element_integpoints*mesh.TotalVolumes;
        }break;
    }

    cout << "Computing interpolants..." << endl;
    vector<KERNEL> kernels;
    switch(mesh.dim)
    {
        case 2:
        {kernels = computeInterpolants(mesh.dim, mesh.points, mesh.surfaces, mesh.TotalSurfaces, mesh.SurfaceNds, gauss);}break;
        case 3:
        {kernels = computeInterpolants(mesh.dim, mesh.points, mesh.volumes, mesh.TotalVolumes, mesh.VolumeNds, gauss);}break;
    }

    cout << "Starting Transient Solver..." << endl;
    int cntr = 1;
    for (int t = 1; t <= settings.nt; t++)
    {
        cout << "Processing Time Step " << t << " of " << settings.nt <<  endl;

        cout << "Assembling global matrices..." << endl;
        GLOBAL_MATRICES global_matrices;
        switch(mesh.dim)
        {
            case 2:
		    {
                global_matrices.assembleDiffusionMatrix(settings.SUPG, mesh.TotalPoints, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts, materials.D, U, zero, zero, settings.stab_prmtr);
                global_matrices.assembleMassMatrix(settings.SUPG, mesh.TotalPoints, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts, U, zero, zero, settings.stab_prmtr);
                global_matrices.assembleForceMatrix(settings.SUPG, mesh.TotalPoints, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts, U, zero, zero, zero, settings.stab_prmtr);
		    }break;
		    case 3:
		    {
                global_matrices.assembleDiffusionMatrix(settings.SUPG, mesh.TotalPoints, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts, materials.D, U, zero, zero, settings.stab_prmtr);
                global_matrices.assembleMassMatrix(settings.SUPG, mesh.TotalPoints, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts, U, zero, zero, settings.stab_prmtr);
                global_matrices.assembleForceMatrix(settings.SUPG, mesh.TotalPoints, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts, U, zero, zero, zero, settings.stab_prmtr);
		    }break;
        }
        cout << "Applying boundary conditions..." << endl;
        global_matrices.applyDBCDirectTransient(mesh.TotalPoints, BC.U.points, BC.U.values, BC.U.total);

        cout << "Solving the system of equations..." << endl;
        U = advanceInTime(settings.transient, U, global_matrices, mesh.TotalPoints, settings.dt, settings.CGiter, settings.tol);

        if (t/settings.prnt_freq == cntr)
        {
            cntr++;
            cout << "Printing Solutions to File..." << endl;
            mesh.printMeshVTK(U, "U", t);
			mesh.printTXT(U,"U",t);
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}

void solveBurgerEBE(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> U)
{
    cout << "SOLVER START..." << endl;
    pltctrl(mesh.points, mesh.TotalPoints, settings.nt, settings.prnt_freq);
    if (mesh.dim == 2)
    {printElements(mesh.surfaces, mesh.SurfaceNds, mesh.TotalSurfaces, "elements");}
    vector<double> zero(mesh.TotalPoints,0);   
    
	cout << "Collect element support for each node..." << endl;
    vector<NBR> ndkernels = computeNodeKernels(mesh.TotalPoints, mesh.volumes, mesh.TotalVolumes, mesh.VolumeNds);

    cout << "Computing quadrature points..." << endl;
    GAUSS gauss;
    int TotQuadPnts = 0;
    switch(mesh.dim)
    {
        case 2:
        {
            gauss.generatePnts(mesh.dim, mesh.SurfaceShape, mesh.ElementOrder);
            TotQuadPnts = gauss.element_integpoints*mesh.TotalSurfaces;
        }break;
        case 3:
        {
            gauss.generatePnts(mesh.dim, mesh.VolumeShape, mesh.ElementOrder);
            TotQuadPnts = gauss.element_integpoints*mesh.TotalVolumes;
        }break;
    }

    cout << "Computing interpolants..." << endl;
    vector<KERNEL> kernels;
    switch(mesh.dim)
    {
        case 2:
        {kernels = computeInterpolants(mesh.dim, mesh.points, mesh.surfaces, mesh.TotalSurfaces, mesh.SurfaceNds, gauss);}break;
        case 3:
        {kernels = computeInterpolants(mesh.dim, mesh.points, mesh.volumes, mesh.TotalVolumes, mesh.VolumeNds, gauss);}break;
    }

    cout << "Collecting all boundary nodes and gauss points with boundary conditions..." << endl;
    vector<int>g_bc; int tot_g_bc=0; vector<int> g_bc2; int tot_g_bc2=0;
    switch(mesh.dim)
    {
        case 2:
        {tie(g_bc, tot_g_bc, g_bc2, tot_g_bc2) = collectGaussPointsWithBC(BC.U.points, BC.U.total, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts);}break;
        case 3:
        {tie(g_bc, tot_g_bc, g_bc2, tot_g_bc2) = collectGaussPointsWithBC(BC.U.points, BC.U.total, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts);}break;
    }

    cout << "Starting Transient Solver..." << endl;
    int cntr = 1;
    for (int t = 1; t <= settings.nt; t++)
    {
        cout << "Processing Time Step " << t << " of " << settings.nt <<  endl;
        switch(mesh.dim)
        {
            case 2:
            {U = solveElementByElement(U, U, zero, zero, zero, materials.D, mesh.TotalPoints, kernels, TotQuadPnts, mesh.surfaces, mesh.SurfaceNds, ndkernels, settings.dt, settings.transient, settings.SUPG, settings.stab_prmtr, settings.CGiter, settings.tol, BC.U.points, BC.U.values, BC.U.total, g_bc, tot_g_bc, g_bc2, tot_g_bc2);}break;
            case 3:
            {U = solveElementByElement(U, U, zero, zero, zero, materials.D, mesh.TotalPoints, kernels, TotQuadPnts, mesh.volumes, mesh.VolumeNds, ndkernels, settings.dt, settings.transient, settings.SUPG, settings.stab_prmtr, settings.CGiter, settings.tol, BC.U.points, BC.U.values, BC.U.total, g_bc, tot_g_bc, g_bc2, tot_g_bc2);}break;
        }
        
        U = setVector(U, BC.U);
        if (t/settings.prnt_freq == cntr)
        {
            cntr++;
            mesh.printMeshVTK(U, "U", t);
			mesh.printTXT(U,"U",t);
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}
