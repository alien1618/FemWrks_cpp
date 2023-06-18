#include "../femwrks.h"

void solveNavierStokes(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Qx, vector<double> Qy, vector<double> Qz)
{
    cout << "SOLVER START..." << endl;
    pltctrl(mesh.points, mesh.TotalPoints, settings.nt, settings.prnt_freq);
    if (mesh.dim == 2)
    {printElements(mesh.surfaces, mesh.SurfaceNds, mesh.TotalSurfaces, "elems");}
    
    cout << "Computing quadrature points..." << endl;
    int TotQuadPnts = 0;
    GAUSS gauss;
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
    vector<double> V(mesh.TotalPoints+1);
    vector<double> P(mesh.TotalPoints+1);
    vector<double> dPdx(mesh.TotalPoints+1);
    vector<double> dPdy(mesh.TotalPoints+1);
    vector<double> dPdz(mesh.TotalPoints+1);
    for (int t = 1; t <= settings.nt; t++)
    {
        cout << "Processing Time Step " << t << " of " << settings.nt <<  endl;

		GLOBAL_MATRICES global_matrices;
		switch(mesh.dim)
        {
            case 2:
		    {
				cout << "Assembling global matrices..." << endl;
                global_matrices.assembleDiffusionMatrix(settings.SUPG, mesh.TotalPoints, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts, materials.Nu, Vx, Vy, Vz, settings.stab_prmtr);
                global_matrices.assembleMassMatrix(settings.SUPG, mesh.TotalPoints, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts, Vx, Vy, Vz, settings.stab_prmtr);

				cout << "Calculating Intermediate Vx..." << endl;
                global_matrices.F = Qx;
				global_matrices.applyDBCDirectTransient(mesh.TotalPoints, BC.Vx.points, BC.Vx.values, BC.Vx.total);
				Vx = advanceInTime(settings.transient, Vx, global_matrices, mesh.TotalPoints, settings.dt, settings.CGiter, settings.tol);

				cout << "Calculating Intermediate Vy..." << endl;
				global_matrices.F = Qy;
				global_matrices.applyDBCDirectTransient(mesh.TotalPoints, BC.Vy.points, BC.Vy.values, BC.Vy.total);
				Vy = advanceInTime(settings.transient, Vy, global_matrices, mesh.TotalPoints, settings.dt, settings.CGiter, settings.tol);

				cout << "Solving the Poisson Equation for Pressure..." << endl;
				tie(P, dPdx, dPdy, dPdz) = computePressure(BC.P, settings.dt, Vx, Vy, Vz,  mesh.TotalPoints, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts, materials.ro, settings.CGiter, settings.tol);
		    }break;
		    case 3:
		    {
				cout << "Assembling global matrices..." << endl;
                global_matrices.assembleDiffusionMatrix(settings.SUPG, mesh.TotalPoints, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts, materials.Nu, Vx, Vy, Vz, settings.stab_prmtr);
                global_matrices.assembleMassMatrix(settings.SUPG, mesh.TotalPoints, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts, Vx, Vy, Vz, settings.stab_prmtr);

				cout << "Calculating Intermediate Vx..." << endl;
                global_matrices.F = Qx;
				global_matrices.applyDBCDirectTransient(mesh.TotalPoints, BC.Vx.points, BC.Vx.values, BC.Vx.total);
				Vx = advanceInTime(settings.transient, Vx, global_matrices, mesh.TotalPoints, settings.dt, settings.CGiter, settings.tol);

				cout << "Calculating Intermediate Vy..." << endl;
                global_matrices.F = Qy;
				global_matrices.applyDBCDirectTransient(mesh.TotalPoints, BC.Vy.points, BC.Vy.values, BC.Vy.total);
				Vy = advanceInTime(settings.transient, Vy, global_matrices, mesh.TotalPoints, settings.dt, settings.CGiter, settings.tol);

				cout << "Calculating Intermediate Vz..." << endl;
                global_matrices.F = Qz;
				global_matrices.applyDBCDirectTransient(mesh.TotalPoints,  BC.Vz.points, BC.Vz.values, BC.Vz.total);
				Vz = advanceInTime(settings.transient, Vz, global_matrices, mesh.TotalPoints, settings.dt, settings.CGiter, settings.tol);

				cout << "Solving the Poisson Equation for Pressure..." << endl;
				tie(P, dPdx, dPdy, dPdz) = computePressure(BC.P, settings.dt, Vx, Vy, Vz, mesh.TotalPoints, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts, materials.ro, settings.CGiter, settings.tol);
		    }break;
        }
		cout << "Computing Corrected Flow Velocities Vx, Vy, Vz..." << endl;
		//#pragma omp parallel for
		for(int i = 1; i <= mesh.TotalPoints; i++)
		{
			Vx[i] = Vx[i]-dPdx[i]*(settings.dt/materials.ro);
			Vy[i] = Vy[i]-dPdy[i]*(settings.dt/materials.ro);
			Vz[i] = Vz[i]-dPdz[i]*(settings.dt/materials.ro);
			V[i] = pow((Vx[i] * Vx[i] + Vy[i] * Vy[i] + Vz[i] * Vz[i]),0.5);
		}
		Vx = setVector(Vx, BC.Vx);
		Vy = setVector(Vy, BC.Vy);
		Vz = setVector(Vz, BC.Vz);
      
        if (t/settings.prnt_freq == cntr)
        {
            cntr++;
            cout << "Printing Solutions to File..." << endl;
            mesh.printMeshVTK(V, "U", t);
            mesh.printMeshVTK(P, "P", t);
            mesh.printVectorsVTK(Vx,Vy,Vz, "vectors", t);
			mesh.printTXT(V,"U",t);
			mesh.printTXT(P,"P",t);
        }
    }
    cout << "SOLVER COMPLETE" << endl;
}

tuple< vector<double>, vector<double>, vector<double>, vector<double> > computePressure(BOUNDARY_CONDITION BC_P, double dt, vector<double> Vx_inter, vector<double> Vy_inter, vector<double> Vz_inter, int TotalPoints, vector<vector<int>> elements, int elem_nds, vector<KERNEL> interpolants, int gauss_TotalPoints, double density, int CGiter, double tol)
{
    GLOBAL_MATRICES G;
    G.K.resize(TotalPoints+1, vector<double> (TotalPoints+1));
    G.F.resize(TotalPoints+1);
    int cmsize = elem_nds+1;
    
    //cout << "looping over all gauss points..." << endl;
//#pragma omp parallel for
    for (int g = 1; g <= gauss_TotalPoints; g++)
    {
        int e = interpolants[g].elem_num;
       
        //cout << "compute gauss matrices" << endl;
        vector<vector<double> > K_gauss(cmsize,vector<double> (cmsize));
        vector<double>  F_gauss(cmsize);

        //donot change...pressure will not be calculated correctly if you modify.
        double dVxdx = 0;
        double dVydy = 0;
        double dVzdz = 0;
        for (int j = 1 ; j <= elem_nds; j++)
        {
            dVxdx = dVxdx+interpolants[g].dNdx[j]*Vx_inter[elements[e][j]];
            dVydy = dVydy+interpolants[g].dNdy[j]*Vy_inter[elements[e][j]];
            dVzdz = dVzdz+interpolants[g].dNdz[j]*Vz_inter[elements[e][j]];
        }

        for (int j = 1 ; j <= elem_nds; j++)
        {
            for (int i = 1; i <= elem_nds; i++)
            {K_gauss[i][j] = (1/density)*(interpolants[g].dNdNx[i][j] + interpolants[g].dNdNy[i][j] + interpolants[g].dNdNz[i][j]);}
            F_gauss[j] = (-1/dt)*(dVxdx + dVydy + dVzdz);
        }

        //cout << "ASSEMBLE THE GAUSS MATRICES DIRECTLY INTO THE GLOBAL MATRICES" << endl;
        for (int n = 1; n <= elem_nds; n++)
        {   int cn = elements[e][n];
            for (int m = 1; m <= elem_nds; m++)
            {
                int cm = elements[e][m];
                G.K[cm][cn] = G.K[cm][cn] + interpolants[g].quad_weight * (K_gauss[m][n]);
            }
            G.F[cn] = G.F[cn] + interpolants[g].quad_weight * F_gauss[n];
        }
    }
    G.applyDBCDirectSteadyState(1, TotalPoints, BC_P.points, BC_P.values, BC_P.total);

	vector<double> P(TotalPoints+1);
	vector<double> dPdx(TotalPoints+1);
 	vector<double> dPdy(TotalPoints+1);   
	vector<double> dPdz(TotalPoints+1);
	P = solvePCG(G.K, G.F, TotalPoints, CGiter, tol);

    // cout << "looping over all gauss points..." << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {
        dPdx[i] = 0;
        dPdy[i] = 0;
        dPdz[i] = 0;
    }
        
//#pragma omp parallel for
    for (int i = 1; i <= gauss_TotalPoints; i++)
    {
        //cout << "compute gauss matrices" << endl;
         int e = interpolants[i].elem_num;
        //donot change...pressure will not be calculated correctly if you modify.
        double dpdx = 0;
        double dpdy = 0;
        double dpdz = 0;
        for (int j = 1 ; j <= elem_nds; j++)
        {
            int nbr = elements[e][j];
            dpdx = dpdx+interpolants[i].dNdx[j]*(P[nbr]);
            dpdy = dpdy+interpolants[i].dNdy[j]*(P[nbr]);
            dpdz = dpdz+interpolants[i].dNdz[j]*(P[nbr]);
        }

        //cout << "ASSEMBLE THE GAUSS MATRICES DIRECTLY INTO THE GLOBAL MATRICES" << endl;
        for (int n = 1; n <= elem_nds; n++)
        {
            int cn = elements[e][n];
            dPdx[cn] = dPdx[cn] + interpolants[i].quad_weight * dpdx;
            dPdy[cn] = dPdy[cn] + interpolants[i].quad_weight * dpdy;
            dPdz[cn] = dPdz[cn] + interpolants[i].quad_weight * dpdz;
        }
    }
	return make_tuple(P, dPdx, dPdy, dPdz);
}
void solveNavierStokesEBE(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Fx, vector<double> Fy, vector<double> Fz)
{
    cout << "SOLVER START..." << endl;
    pltctrl(mesh.points, mesh.TotalPoints, settings.nt, settings.prnt_freq);
    if (mesh.dim == 2)
    {printElements(mesh.surfaces, mesh.SurfaceNds, mesh.TotalSurfaces, "elems");}
    
    cout << "Initializing Variables..." << endl;
    Vx = setVector(Vx, BC.Vx);
    Vy = setVector(Vy, BC.Vy);
    Vz = setVector(Vz, BC.Vz);
    vector<double> V(mesh.TotalPoints+1);
    int cntr = 1;
    
    cout << "Collect element support for each node..." << endl;
    vector<NBR> ndkernels = computeNodeKernels(mesh.TotalPoints, mesh.volumes, mesh.TotalVolumes, mesh.VolumeNds);

    cout << "Computing quadrature points..." << endl;
    int TotQuadPnts = 0;
    GAUSS gauss;
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
    vector<int>g_bc_vx; int tot_g_bc_vx=0; vector<int> g_bc2_vx; int tot_g_bc2_vx=0;
    vector<int>g_bc_vy; int tot_g_bc_vy=0; vector<int> g_bc2_vy; int tot_g_bc2_vy=0;
    vector<int>g_bc_vz; int tot_g_bc_vz=0; vector<int> g_bc2_vz; int tot_g_bc2_vz=0;
    vector<int>g_bc_p; int tot_g_bc_p=0; vector<int> g_bc2_p; int tot_g_bc2_p=0;
    switch(mesh.dim)
    {
        case 2:
    	{
			tie(g_bc_vx, tot_g_bc_vx, g_bc2_vx, tot_g_bc2_vx) = collectGaussPointsWithBC(BC.Vx.points, BC.Vx.total, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts);
			tie(g_bc_vy, tot_g_bc_vy, g_bc2_vy, tot_g_bc2_vy) = collectGaussPointsWithBC(BC.Vy.points, BC.Vy.total, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts);
			tie(g_bc_vz, tot_g_bc_vz, g_bc2_vz, tot_g_bc2_vz) = collectGaussPointsWithBC(BC.Vz.points, BC.Vz.total, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts);
			tie(g_bc_p, tot_g_bc_p, g_bc2_p, tot_g_bc2_p) = collectGaussPointsWithBC(BC.P.points, BC.P.total, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts);
		}break;
        case 3:
    	{
			tie(g_bc_vx, tot_g_bc_vx, g_bc2_vx, tot_g_bc2_vx) = collectGaussPointsWithBC(BC.Vx.points, BC.Vx.total, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts);
			tie(g_bc_vy, tot_g_bc_vy, g_bc2_vy, tot_g_bc2_vy) = collectGaussPointsWithBC(BC.Vy.points, BC.Vy.total, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts);
			tie(g_bc_vz, tot_g_bc_vz, g_bc2_vz, tot_g_bc2_vz) = collectGaussPointsWithBC(BC.Vz.points, BC.Vz.total, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts);
			tie(g_bc_p, tot_g_bc_p, g_bc2_p, tot_g_bc2_p) = collectGaussPointsWithBC(BC.P.points, BC.P.total, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts);
		}break;
    }
    
	cout << "Starting Transient Solver..." << endl;
    vector<double> Vx_iter(mesh.TotalPoints+1);
    vector<double> Vy_iter(mesh.TotalPoints+1);
    vector<double> Vz_iter(mesh.TotalPoints+1);
    vector<double> V_iter(mesh.TotalPoints+1);
    vector<double> V_iter_old(mesh.TotalPoints+1);
    vector<double> P(mesh.TotalPoints+1);
    vector<double> dPdx(mesh.TotalPoints+1);
    vector<double> dPdy(mesh.TotalPoints+1);
    vector<double> dPdz(mesh.TotalPoints+1);
    Vx_iter = Vx;
    Vy_iter = Vy;
    Vz_iter = Vz;
    V_iter = V;
    for (int t = 1; t <= settings.nt; t++)
    {
        cout << "Processing Time Step " << t << " of " << settings.nt << " Complete" << endl;
        for (int s = 1; s <= 1; s++)
        {
			switch(mesh.dim)
            {
                case 2:
				{
					Vx_iter = solveElementByElement(Vx, Vx_iter, Vy_iter, Vz_iter, Fx, materials.Nu, mesh.TotalPoints, kernels, TotQuadPnts, mesh.surfaces, mesh.SurfaceNds, ndkernels, settings.dt, settings.transient, settings.SUPG, settings.stab_prmtr, settings.CGiter, settings.tol, BC.Vx.points, BC.Vx.values, BC.Vx.total, g_bc_vx, tot_g_bc_vx, g_bc2_vx, tot_g_bc2_vx);
			
					Vy_iter = solveElementByElement(Vy, Vx_iter, Vy_iter, Vz_iter, Fy, materials.Nu, mesh.TotalPoints, kernels, TotQuadPnts, mesh.surfaces, mesh.SurfaceNds, ndkernels, settings.dt, settings.transient, settings.SUPG, settings.stab_prmtr, settings.CGiter, settings.tol, BC.Vy.points, BC.Vy.values, BC.Vy.total, g_bc_vy, tot_g_bc_vy, g_bc2_vy, tot_g_bc2_vy);
				
					tie(P, dPdx, dPdy, dPdz) = computePressureEBE(Vx_iter, Vy_iter, Vz_iter, mesh.TotalPoints, kernels, TotQuadPnts,mesh.surfaces, mesh.SurfaceNds, ndkernels, materials.ro, settings.dt, settings.CGiter, settings.tol, BC.P.points, BC.P.values, BC.P.total, g_bc_p, tot_g_bc_p, g_bc2_p, tot_g_bc2_p);
				}break;
                case 3:
            	{
					Vx_iter = solveElementByElement(Vx, Vx_iter, Vy_iter, Vz_iter, Fx, materials.Nu, mesh.TotalPoints, kernels, TotQuadPnts, mesh.volumes, mesh.VolumeNds, ndkernels, settings.dt, settings.transient, settings.SUPG, settings.stab_prmtr, settings.CGiter, settings.tol, BC.Vx.points, BC.Vx.values, BC.Vx.total, g_bc_vx, tot_g_bc_vx, g_bc2_vx, tot_g_bc2_vx);
			
					Vy_iter = solveElementByElement(Vy, Vx_iter, Vy_iter, Vz_iter, Fy, materials.Nu, mesh.TotalPoints, kernels, TotQuadPnts, mesh.volumes, mesh.VolumeNds, ndkernels, settings.dt, settings.transient, settings.SUPG, settings.stab_prmtr, settings.CGiter, settings.tol, BC.Vy.points, BC.Vy.values, BC.Vy.total, g_bc_vy, tot_g_bc_vy, g_bc2_vy, tot_g_bc2_vy);
				
					Vz_iter = solveElementByElement(Vz, Vx_iter, Vy_iter,Vz_iter, Fz, materials.Nu, mesh.TotalPoints, kernels, TotQuadPnts, mesh.volumes, mesh.VolumeNds, ndkernels, settings.dt, settings.transient, settings.SUPG, settings.stab_prmtr, settings.CGiter, settings.tol, BC.Vz.points, BC.Vz.values, BC.Vz.total, g_bc_vz, tot_g_bc_vz, g_bc2_vz, tot_g_bc2_vz);
			    
			    	tie(P, dPdx, dPdy, dPdz) = computePressureEBE(Vx_iter, Vy_iter, Vz_iter, mesh.TotalPoints, kernels, TotQuadPnts,mesh.volumes, mesh.VolumeNds, ndkernels, materials.ro, settings.dt, settings.CGiter, settings.tol, BC.P.points, BC.P.values, BC.P.total, g_bc_p, tot_g_bc_p, g_bc2_p, tot_g_bc2_p);
				}break;
            }
               
            //cout << "Computing Corrected Flow Velocities Vx, Vy..." << endl;
            //#pragma omp parallel for
            V_iter_old = V_iter;
            for(int i = 1; i <= mesh.TotalPoints; i++)
            {
                Vx_iter[i] = Vx_iter[i]-dPdx[i]*(settings.dt/materials.ro);
                Vy_iter[i] = Vy_iter[i]-dPdy[i]*(settings.dt/materials.ro);
                Vz_iter[i] = Vz_iter[i]-dPdz[i]*(settings.dt/materials.ro);
            }
            Vx_iter = setVector(Vx_iter, BC.Vx);
            Vy_iter = setVector(Vy_iter, BC.Vy);
            Vz_iter = setVector(Vz_iter, BC.Vz);
            double res = 0;
            for (int i = 1; i <= mesh.TotalPoints; i++)
            {
                V_iter[i] = pow((Vx_iter[i] * Vx_iter[i] + Vy_iter[i] * Vy_iter[i]+ Vz_iter[i] * Vz_iter[i]),0.5);
                res = res+abs(V_iter[i]-V_iter_old[i]);
            }
            cout << "iteration " << s << " RES = " << res << endl;
            if (t/settings.prnt_freq == cntr)
            {
                cntr++;
				mesh.printMeshVTK(V, "U", t);
				mesh.printMeshVTK(P, "P", t);
				mesh.printVectorsVTK(Vx,Vy,Vz, "vectors", t);
				mesh.printTXT(V,"U",t);
				mesh.printTXT(P,"P",t);
            }
        }
        V = V_iter;
        Vx = Vx_iter;
        Vy = Vy_iter;
    }
    cout << "SOLVER COMPLETE" << endl;
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>> computePressureEBE(vector<double> Vx_inter, vector<double> Vy_inter, vector<double> Vz_inter, int TotalMeshPoints, vector<KERNEL> kernels, int gauss_TotalPoints, vector<vector<int>> elements, int ElementNds, vector<NBR> ndkernels, double ro, double dt, int CGiter, double conv_tol, vector<int> bc_ndnum, vector<double> bc_ndval, int BC_TotNds, vector<int> g_bc, int tot_g_bc, vector<int> g_bc2, int tot_g_bc2)
{
    vector<double> P(TotalMeshPoints+1);
    vector<double> dPdx(TotalMeshPoints+1);
    vector<double> dPdy(TotalMeshPoints+1);
    vector<double> dPdz(TotalMeshPoints+1);
    vector<double> F(TotalMeshPoints+1);

    //cout << "looping over all gauss points..." << endl;
    for (int g = 1; g <= gauss_TotalPoints; g++)
    {
        int e = kernels[g].elem_num;
        
		//cout << "compute gauss matrices" << endl;
        vector<vector<double> > K_gauss(ElementNds+1,vector<double> (ElementNds+1));

        //donot change...pressure will not be calculated correctly if you modify.
        double dVxdx = 0;
        double dVydy = 0;
        double dVzdz = 0;
        for (int j = 1 ; j <= ElementNds; j++)
        {
            dVxdx = dVxdx+kernels[g].dNdx[j]*Vx_inter[elements[e][j]];
            dVydy = dVydy+kernels[g].dNdy[j]*Vy_inter[elements[e][j]];
            dVzdz = dVzdz+kernels[g].dNdz[j]*Vz_inter[elements[e][j]];
        }

        for (int j = 1 ; j <= ElementNds; j++)
        {
            for (int i = 1; i <= ElementNds; i++)
            {K_gauss[i][j] = kernels[g].quad_weight*(kernels[g].dNdNx[i][j] + kernels[g].dNdNy[i][j] + kernels[g].dNdNz[i][j]);}
            int cj = elements[e][j];
            F[cj] = F[cj] + kernels[g].quad_weight*(-ro/dt)*(dVxdx + dVydy + dVzdz);
        }
        kernels[g].MK = K_gauss;
    }
    P = solvePCG(1, F, ndkernels, TotalMeshPoints, elements, ElementNds, kernels, bc_ndnum, bc_ndval, BC_TotNds, g_bc, tot_g_bc, g_bc2, tot_g_bc2, CGiter, conv_tol);

    // cout << "computing the pressure gradients..." << endl;
    for (int i = 1; i <= TotalMeshPoints; i++)
    {
        dPdx[i] = 0;
        dPdy[i] = 0;
        dPdz[i] = 0;
    }
    for (int g = 1; g <= gauss_TotalPoints; g++)
    {
        //cout << "compute gauss matrices" << endl;
         int e = kernels[g].elem_num;
		//donot change...pressure will not be calculated correctly if you modify.
        double dpdx = 0;
        double dpdy = 0;
        double dpdz = 0;
        for (int j = 1 ; j <= ElementNds; j++)
        {
            int nbr = elements[e][j];
            dpdx = dpdx+kernels[g].dNdx[j]*(P[nbr]);
            dpdy = dpdy+kernels[g].dNdy[j]*(P[nbr]);
            dpdz = dpdz+kernels[g].dNdz[j]*(P[nbr]);
        }

        //cout << "ASSEMBLE THE GAUSS MATRICES DIRECTLY INTO THE GLOBAL MATRICES" << endl;
        for (int n = 1; n <= ElementNds; n++)
        {
            int cn = elements[e][n];
            dPdx[cn] = dPdx[cn] + kernels[g].quad_weight * dpdx;
            dPdy[cn] = dPdy[cn] + kernels[g].quad_weight * dpdy;
            dPdz[cn] = dPdz[cn] + kernels[g].quad_weight * dpdz;
        }
    }
    return make_tuple(P, dPdx, dPdy,dPdz);
}
