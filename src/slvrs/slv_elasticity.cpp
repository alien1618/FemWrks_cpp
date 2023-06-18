#include "../femwrks.h"

vector<vector<double> > computeD2D(int load_type, double E, double nu, double t)
{
	vector<vector<double> > D_matrix(4, vector<double>(4));
	double D1 = 1;
	double D2 = 2;
	double D3 = 3;
	switch(load_type)
    {
        case 1: // plane stress
	    {
		    D1 = E*t/(1-(nu*nu));
		    D2 = nu*D1;
		    D3 = D1*(1-nu)/2;
	    }break;
        case 2: // plane strain
	    {
		    D1 = (1-nu)*E*t/((1+nu)*(1-2*nu));
		    D2 = nu*E*t/((1+nu)*(1-2*nu));
		    D3 = (1-2*nu)*E*t/2*((1+nu)*(1-2*nu));
	    }break;
    }
	D_matrix[1][1] = D1; D_matrix[2][1] = D2;D_matrix[3][1] = 0;
	D_matrix[1][2] = D2; D_matrix[2][2] = D1; D_matrix[3][2] = 0;
	D_matrix[1][3] = 0; D_matrix[2][3] = 0; D_matrix[3][3] = D3;

	return D_matrix;
}
vector<vector<double> > computeD3D(double E, double nu)
{
	vector<vector<double> > D_matrix(7, vector<double>(7));
	double t1 = (E/((1+nu)*(1-2*nu)));
	double D1 = (1-nu)*t1;
	double D2 = nu*t1;
	double D3 = 0.5*(1-2*nu)*t1;
	D_matrix[1][1] = D1; D_matrix[2][1] = D2; D_matrix[3][1] = D2;
	D_matrix[1][2] = D2; D_matrix[2][2] = D1; D_matrix[3][2] = D2;
	D_matrix[1][3] = D2; D_matrix[2][3] = D2; D_matrix[3][3] = D1;
	D_matrix[4][4] = D3; D_matrix[5][5] = D3; D_matrix[6][6] = D3;

	return D_matrix;
}
vector<vector<double> > computeB3D(int VolumeNds, vector<double> dNdx, vector<double> dNdy, vector<double> dNdz)
{
	int k = 1;
	int DOF = 3;
	vector<vector<double> > B(DOF*VolumeNds+1, vector<double>(7));
	for (int i = 1; i <= VolumeNds; i++)
	{
		B[k][1] = dNdx[i];
		B[DOF*i-1][2] = dNdy[i];
		B[DOF*i][3] = dNdz[i];
		B[k][4] = dNdy[i];
		B[k+1][4] = dNdx[i];
		B[k+1][5] = dNdz[i];
		B[k+2][5] = dNdy[i];
		B[k][6] = dNdz[i];
		B[k+2][6] = dNdx[i];
		k = k + 3;
	}
	return B;
}
vector<vector<double> > computeB2D(int SurfaceNds, vector<double> dNdx, vector<double> dNdy)
{
	int k = 1;
	int DOF = 2;
	vector<vector<double> > B(DOF*SurfaceNds+1, vector<double>(4));
	for (int i = 1; i <= SurfaceNds; i++)
	{
		B[k][1] = dNdx[i];
		B[DOF*i][2] = dNdy[i];
		B[DOF*i][3] = dNdx[i];
		B[k][3] = dNdy[i];
		k = k + 2;
	}
	return B;
}
void solveElasticity(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings)
{
    cout << "SOLVER START..." << endl;
    pltctrl(mesh.points, mesh.TotalPoints, settings.nt, settings.prnt_freq);
    if (mesh.dim == 2)
    {printElements(mesh.surfaces, mesh.SurfaceNds, mesh.TotalSurfaces, "elems");}
    
    cout << "Initializing variables..." << endl;
    int DOF = mesh.dim;
    vector<double> U(mesh.TotalPoints*DOF);
    mesh.printMeshVTK(U, "U", 0);
	mesh.printTXT(U, "U", 0);
    
	cout << "Collecting support elements around each node..." << endl;   
	vector<NBR> ndkernels;
    switch(mesh.dim)
    {
        case 2:
	    {ndkernels = computeNodeKernels(mesh.TotalPoints, mesh.surfaces, mesh.TotalSurfaces, mesh.SurfaceNds);}break;
        case 3:
	    {ndkernels = computeNodeKernels(mesh.TotalPoints, mesh.volumes, mesh.TotalVolumes, mesh.VolumeNds);}break;
    }
	
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

    cout << "Computing local matrices and assembling global matrices..." << endl;
    GLOBAL_MATRICES global_matrices;
    switch(mesh.dim)
    {
        case 2:
		{
            global_matrices.assembleStiffnessMatrix(mesh.dim, 1, materials.MOE, materials.nu, materials.thickness, mesh.TotalPoints, kernels, TotQuadPnts,  mesh.surfaces, mesh.SurfaceNds);
			cout << "Applying boundary conditions..." << endl;
			int iter = 1;
			for (int j = 1; j <= BC.Fx.total; j++)
			{
				int ndnum = BC.Fx.points[j];
				global_matrices.F[DOF*ndnum-iter] = BC.Fx.values[j];
			}
			for (int j = 1; j <= BC.Fy.total; j++)
			{
				int ndnum = BC.Fy.points[j];
				global_matrices.F[DOF*ndnum] = BC.Fy.values[j];
			}
		}break;
        case 3:
		{
            global_matrices.assembleStiffnessMatrix(mesh.dim, 1, materials.MOE, materials.nu, materials.thickness, mesh.TotalPoints, kernels, TotQuadPnts,  mesh.volumes, mesh.VolumeNds);
			cout << "Applying boundary conditions..." << endl;
			int iter = 2;
			for (int j = 1; j <= BC.Fx.total; j++)
			{
				int ndnum = BC.Fx.points[j];
				global_matrices.F[DOF*ndnum-iter] = BC.Fx.values[j];
			}
			iter = 1;
			for (int j = 1; j <= BC.Fy.total; j++)
			{
				int ndnum = BC.Fy.points[j];
				global_matrices.F[DOF*ndnum-1] = BC.Fy.values[j];
			}
			for (int j = 1; j <= BC.Fz.total; j++)
			{
				int ndnum = BC.Fz.points[j];
				global_matrices.F[DOF*ndnum] = BC.Fz.values[j];
			}
		}break;
    }
    global_matrices.applyDBCDirectSteadyState(DOF, mesh.TotalPoints, BC.U.points, BC.U.values, BC.U.total);

    cout << "Solving the system of equations..." << endl;
    U = solvePCG(global_matrices.K,  global_matrices.F, DOF*mesh.TotalPoints, 1000, 0.0001);
	
	cout << "Computing stress distribution..." << endl;
	vector<double> deformation(mesh.TotalPoints+1);
	vector<double> vms(mesh.TotalPoints+1);
	vector<double> FOS(mesh.TotalPoints+1);
	vector<double> gauss_vms(TotQuadPnts+1);	
	
    //collecting all quadrature points inside each element
    vector<NBR> elems;
    switch(mesh.dim)
    {
        case 3:
	    {
		    elems.resize(mesh.TotalVolumes+1);
		    for (int i = 1; i <= mesh.TotalVolumes; i++)
		    {
			    elems[i].TotalNbrs = 0;
			    elems[i].nbrs.push_back(1);
		    }
	    }break;
        case 2:
        {
		    elems.resize(mesh.TotalSurfaces+1);
		    for (int i = 1; i <= mesh.TotalSurfaces; i++)
		    {
			    elems[i].TotalNbrs = 0;
			    elems[i].nbrs.push_back(1);
		    }
	    }break;
    }
	for (int i = 1; i <= TotQuadPnts; i++)
	{
		elems[kernels[i].elem_num].TotalNbrs++;
		elems[kernels[i].elem_num].nbrs.push_back(1);
		elems[kernels[i].elem_num].nbrs[elems[kernels[i].elem_num].TotalNbrs] = i;
	}
	switch(mesh.dim)
    {
        case 3:
	    {
			vector<double> U_supp(DOF*mesh.VolumeNds+1);
			vector<double> dof(DOF*mesh.VolumeNds+1);
			vector<vector<double> > D_matrix = computeD3D(materials.MOE, materials.nu);
			for (int g = 1; g <= TotQuadPnts; g++)
			{
				int e = kernels[g].elem_num;
				vector<vector<double> > B = computeB3D(mesh.VolumeNds, kernels[g].dNdx, kernels[g].dNdy, kernels[g].dNdz);
				int kk = 1;
				for (int j = 1; j <= mesh.VolumeNds; j++)
				{
					for (int iter = DOF-1; iter >= 0; iter--)
					{
						dof[kk] =  DOF*mesh.volumes[e][j]-iter;
						kk = kk + 1;
					}
				}
				for (int i = 1; i <= DOF*mesh.VolumeNds; i++)
				{U_supp[i] = U[dof[i]];}
		   
				double strain_xx = 0;
				double strain_yy = 0;
				double strain_zz = 0;
				double gamma_xy = 0;
				double gamma_xz = 0;
				double gamma_yz = 0;
				for (int s = 1; s <= DOF*mesh.VolumeNds; s++)
				{
					strain_xx = strain_xx  + B[s][1]*U_supp[s];
					strain_yy = strain_yy  + B[s][2]*U_supp[s];
					strain_zz = strain_zz  + B[s][3]*U_supp[s];
					gamma_xy = gamma_xy  + B[s][4]*U_supp[s];
					gamma_xz = gamma_xz  + B[s][6]*U_supp[s];
					gamma_yz = gamma_yz  + B[s][5]*U_supp[s];
				}
				double stress_xx = D_matrix[1][1]*strain_xx + D_matrix[2][1]*strain_yy + D_matrix[3][1]*strain_zz+D_matrix[4][1]*gamma_xy+D_matrix[5][1]*gamma_xz+D_matrix[6][1]*gamma_yz;
				double stress_yy = D_matrix[1][2]*strain_xx + D_matrix[2][2]*strain_yy + D_matrix[3][2]*strain_zz+D_matrix[4][2]*gamma_xy+D_matrix[5][2]*gamma_xz+D_matrix[6][2]*gamma_yz;
				double stress_zz = D_matrix[1][3]*strain_xx + D_matrix[2][3]*strain_yy + D_matrix[3][3]*strain_zz+D_matrix[4][3]*gamma_xy+D_matrix[5][3]*gamma_xz+D_matrix[6][3]*gamma_yz;

				gauss_vms[g] = sqrt(pow((stress_xx-stress_yy),2)+pow((stress_yy-stress_zz),2)+pow((stress_xx-stress_zz),2)+ 6*(gamma_xy*gamma_xy+gamma_yz*gamma_yz+gamma_xz*gamma_xz))/sqrt(2.0);
				//gauss_vms[g] = sqrt(pow((stress_xx-stress_yy),2)+pow((stress_yy-stress_zz),2)+pow((stress_xx-stress_zz),2)+ 6*(stress_xy*stress_xy+stress_yz+stress_yz+stress_xz*stress_xz))/sqrt(2.0);
			}
			
			cout << "Computing final nodal deformation..." << endl;
			for (int i = 1; i <= mesh.TotalPoints; i++)
			{
				mesh.points[i].x = mesh.points[i].x + U[DOF*i-2];
				mesh.points[i].y = mesh.points[i].y + U[DOF*i-1];
				mesh.points[i].z = mesh.points[i].z + U[DOF*i];
				deformation[i] = sqrt(U[DOF*i-2]*U[DOF*i-2] + U[DOF*i-1]*U[DOF*i-1] + U[DOF*i]*U[DOF*i]);
			}
		}break;
        case 2:
    	{
            int load_type = 1;
            vector<double> U_supp(DOF*mesh.SurfaceNds+1);
            vector<double> dof(DOF*mesh.SurfaceNds+1);
            vector<vector<double> > D_matrix = computeD2D(load_type, materials.MOE, materials.nu, materials.thickness);
            for (int g = 1; g <= TotQuadPnts; g++)
            {
                int e = kernels[g].elem_num;
                vector<vector<double> > B = computeB2D(mesh.SurfaceNds, kernels[g].dNdx, kernels[g].dNdy);
                int kk = 1;
                for (int j = 1; j <= mesh.SurfaceNds; j++)
                {
                    for (int iter = DOF-1; iter >= 0; iter--)
                    {
                        dof[kk] =  DOF*mesh.surfaces[e][j]-iter;
                        kk = kk + 1;
                    }
                }		
                for (int i = 1; i <= DOF*mesh.SurfaceNds; i++)
                {U_supp[i] = U[dof[i]];}
            
                double strain_xx = 0;
                double strain_yy = 0;
                double gamma_xy = 0;
                for (int s = 1; s <= DOF*mesh.SurfaceNds; s++)
                {
                    strain_xx = strain_xx  + B[s][1]*U_supp[s];
                    strain_yy = strain_yy  + B[s][2]*U_supp[s];
                    gamma_xy = gamma_xy  + B[s][3]*U_supp[s];
                }
                double stress_xx = D_matrix[1][1]*strain_xx + D_matrix[2][1]*strain_yy + D_matrix[3][1]*gamma_xy;
                double stress_yy = D_matrix[1][2]*strain_xx + D_matrix[2][2]*strain_yy + D_matrix[3][2]*gamma_xy;
                //stress_xy = D_matrix[1][3]*strain_xx + D_matrix[2][3]*strain_yy + D_matrix[3][3]*gamma_xy;
                gauss_vms[g] = sqrt(pow((stress_xx-stress_yy),2)+ stress_xx*stress_xx + stress_yy*stress_yy + 6*gamma_xy*gamma_xy)/sqrt(2.0);
            }
	
            cout << "Computing final nodal deformation..." << endl;
            for (int i = 1; i <= mesh.TotalPoints; i++)
            {
                mesh.points[i].x = mesh.points[i].x + U[DOF*i-1];
                mesh.points[i].y = mesh.points[i].y + U[DOF*i];
                deformation[i] = sqrt(U[DOF*i-1]*U[DOF*i-1] + U[DOF*i]*U[DOF*i]);
            }
        }break;
    }
	//interpolate stresses from gauss points to nodes
	for (int i = 1; i <= mesh.TotalPoints; i++)
	{
		double sum_w = 0;
		double sum_vms = 0;
		for (int j = 1; j <= ndkernels[i].TotalNbrs; j++)
		{
			int e = ndkernels[i].nbrs[j];
			for (int g = 1; g <= elems[e].TotalNbrs; g++)
			{
				int q = elems[e].nbrs[g];
				double dx = mesh.points[i].x - kernels[q].pnt.x;
				double dy = mesh.points[i].y - kernels[q].pnt.y;
				double dz = mesh.points[i].z - kernels[q].pnt.z;
				double d = sqrt(dx*dx+dy*dy+dz*dz);
				double w = 1/d;
				sum_w = sum_w + w;
				sum_vms = sum_vms + w*gauss_vms[q];
			}
		}
		vms[i] = sum_vms/sum_w;
		FOS[i] = vms[i]/materials.YTS;
	}

	cout << "Printing results to file..." << endl;
	mesh.printMeshVTK(deformation, "U", 1);
	mesh.printMeshVTK(vms,"vms",1);
	mesh.printMeshVTK(FOS,"fos",1);
	
	mesh.printTXT(U, "U", 1);
    mesh.printTXT(vms, "vms", 1);
    mesh.printTXT(FOS, "fos", 1);

    cout << "SOLVER COMPLETE" << endl;   
}
void solveElasticityEBE(MESH mesh, MATERIALS materials, BOUNDARY_CONDITIONS BC, SOLVER_SETTINGS settings)
{
    cout << "SOLVER START..." << endl;
    pltctrl(mesh.points, mesh.TotalPoints, settings.nt, settings.prnt_freq);
    if (mesh.dim == 2)
    {printElements(mesh.surfaces, mesh.SurfaceNds, mesh.TotalSurfaces, "elems");}
    
    cout << "Initializing variables..." << endl;
    int DOF = mesh.dim;
    vector<double> U(mesh.TotalPoints*DOF);
    mesh.printMeshVTK(U, "U", 0);
	mesh.printTXT(U, "U", 0);
    
	cout << "Collecting support elements around each node..." << endl;   
	vector<NBR> ndkernels;
    switch(mesh.dim)
    {
        case 2:
	    {ndkernels = computeNodeKernels(mesh.TotalPoints, mesh.surfaces, mesh.TotalSurfaces, mesh.SurfaceNds);}break;
        case 3:
	    {ndkernels = computeNodeKernels(mesh.TotalPoints, mesh.volumes, mesh.TotalVolumes, mesh.VolumeNds);}break;
    }
    
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
    vector<int>g_bc; int tot_g_bc=0; vector<int> g_bc2; int tot_g_bc2=0;
    switch(mesh.dim)
    {
        case 2:
        {tie(g_bc, tot_g_bc, g_bc2, tot_g_bc2) = collectGaussPointsWithBC(BC.U.points, BC.U.total, mesh.surfaces, mesh.SurfaceNds, kernels, TotQuadPnts);}break;
        case 3:
        {tie(g_bc, tot_g_bc, g_bc2, tot_g_bc2) = collectGaussPointsWithBC(BC.U.points, BC.U.total, mesh.volumes, mesh.VolumeNds, kernels, TotQuadPnts);}break;
    }
	
    vector<double> F(DOF*mesh.TotalPoints+1);
    for(int j = 1; j <= DOF*mesh.TotalPoints; j++)
    {F[j] = 0;}

    cout << "Calculating the force vector..." << endl;
    switch(mesh.dim)
    {
        case 2:
        {
            int iter = 1;
            for (int j = 1; j <= BC.Fx.total; j++)
            {
                int ndnum = BC.Fx.points[j];
                F[DOF*ndnum-iter] = BC.Fx.values[j];
            }
            for (int j = 1; j <= BC.Fy.total; j++)
            {
                int ndnum = BC.Fy.points[j];
                F[DOF*ndnum] = BC.Fy.values[j];
            }
        }break;
        case 3:
        {
            cout << "Applying boundary conditions..." << endl;
            int iter = 2;
            for (int j = 1; j <= BC.Fx.total; j++)
            {
                int ndnum = BC.Fx.points[j];
                F[DOF*ndnum-iter] = BC.Fx.values[j];
            }
            iter = 1;
            for (int j = 1; j <= BC.Fy.total; j++)
            {
                int ndnum = BC.Fy.points[j];
                F[DOF*ndnum-1] = BC.Fy.values[j];
            }
            for (int j = 1; j <= BC.Fz.total; j++)
            {
                int ndnum = BC.Fz.points[j];
                F[DOF*ndnum] = BC.Fz.values[j];
            }
        }break;
    }

	cout << "Solve Element-by-Element..." << endl;
	switch(mesh.dim)
	{
		case 2:
		{
			U = solveElementByElement(mesh.dim, DOF, materials.MOE, materials.nu, materials.thickness, mesh.TotalPoints, kernels, TotQuadPnts, mesh.surfaces, mesh.SurfaceNds, ndkernels, settings.CGiter, settings.tol, F, BC.U.points, BC.U.values, BC.U.total, g_bc, tot_g_bc, g_bc2, tot_g_bc2);
		}break;
		case 3:
		{
			U = solveElementByElement(mesh.dim, DOF, materials.MOE, materials.nu, materials.thickness, mesh.TotalPoints, kernels, TotQuadPnts, mesh.volumes, mesh.VolumeNds, ndkernels, settings.CGiter, settings.tol, F, BC.U.points, BC.U.values, BC.U.total, g_bc, tot_g_bc, g_bc2, tot_g_bc2);
		}break;
	}
	cout << "Computing stress distribution..." << endl;
	vector<double> deformation(mesh.TotalPoints+1);
	vector<double> vms(mesh.TotalPoints+1);
	vector<double> FOS(mesh.TotalPoints+1);
	vector<double> gauss_vms(TotQuadPnts+1);	
	
    //collecting all quadrature points inside each element
    vector<NBR> elems;
    switch(mesh.dim)
    {
        case 3:
	    {
		    elems.resize(mesh.TotalVolumes+1);
		    for (int i = 1; i <= mesh.TotalVolumes; i++)
		    {
			    elems[i].TotalNbrs = 0;
			    elems[i].nbrs.push_back(1);
		    }
	    }break;
        case 2:
        {
		    elems.resize(mesh.TotalSurfaces+1);
		    for (int i = 1; i <= mesh.TotalSurfaces; i++)
		    {
			    elems[i].TotalNbrs = 0;
			    elems[i].nbrs.push_back(1);
		    }
	    }break;
    }
	for (int i = 1; i <= TotQuadPnts; i++)
	{
		elems[kernels[i].elem_num].TotalNbrs++;
		elems[kernels[i].elem_num].nbrs.push_back(1);
		elems[kernels[i].elem_num].nbrs[elems[kernels[i].elem_num].TotalNbrs] = i;
	}
	switch(mesh.dim)
    {
        case 3:
        {
            vector<vector<double> > D_matrix = computeD3D(materials.MOE, materials.nu);
            vector<double> U_supp(DOF*mesh.VolumeNds+1);
            vector<double> dof(DOF*mesh.VolumeNds+1);
            for (int g = 1; g <= TotQuadPnts; g++)
            {
                int e = kernels[g].elem_num;
                vector<vector<double> > B = computeB3D(mesh.VolumeNds, kernels[g].dNdx, kernels[g].dNdy, kernels[g].dNdz);
                int kk = 1;
                for (int j = 1; j <= mesh.VolumeNds; j++)
                {
                    for (int iter = DOF-1; iter >= 0; iter--)
                    {
                        dof[kk] =  DOF*mesh.volumes[e][j]-iter;
                        kk = kk + 1;
                    }
                }
                for (int i = 1; i <= DOF*mesh.VolumeNds; i++)
                {U_supp[i] = U[dof[i]];}
               
                double strain_xx = 0;
                double strain_yy = 0;
                double strain_zz = 0;
                double gamma_xy = 0;
                double gamma_xz = 0;
                double gamma_yz = 0;
                for (int s = 1; s <= DOF*mesh.VolumeNds; s++)
                {
                    strain_xx = strain_xx  + B[s][1]*U_supp[s];
                    strain_yy = strain_yy  + B[s][2]*U_supp[s];
                    strain_zz = strain_zz  + B[s][3]*U_supp[s];
                    gamma_xy = gamma_xy  + B[s][4]*U_supp[s];
                    gamma_xz = gamma_xz  + B[s][6]*U_supp[s];
                    gamma_yz = gamma_yz  + B[s][5]*U_supp[s];
                }
                double stress_xx = D_matrix[1][1]*strain_xx + D_matrix[2][1]*strain_yy + D_matrix[3][1]*strain_zz+D_matrix[4][1]*gamma_xy+D_matrix[5][1]*gamma_xz+D_matrix[6][1]*gamma_yz;
                double stress_yy = D_matrix[1][2]*strain_xx + D_matrix[2][2]*strain_yy + D_matrix[3][2]*strain_zz+D_matrix[4][2]*gamma_xy+D_matrix[5][2]*gamma_xz+D_matrix[6][2]*gamma_yz;
                double stress_zz = D_matrix[1][3]*strain_xx + D_matrix[2][3]*strain_yy + D_matrix[3][3]*strain_zz+D_matrix[4][3]*gamma_xy+D_matrix[5][3]*gamma_xz+D_matrix[6][3]*gamma_yz;

                gauss_vms[g] = sqrt(pow((stress_xx-stress_yy),2)+pow((stress_yy-stress_zz),2)+pow((stress_xx-stress_zz),2)+ 6*(gamma_xy*gamma_xy+gamma_yz*gamma_yz+gamma_xz*gamma_xz))/sqrt(2.0);
                //gauss_vms[g] = sqrt(pow((stress_xx-stress_yy),2)+pow((stress_yy-stress_zz),2)+pow((stress_xx-stress_zz),2)+ 6*(stress_xy*stress_xy+stress_yz+stress_yz+stress_xz*stress_xz))/sqrt(2.0);
            }
                
            cout << "Computing final nodal deformation..." << endl;
            for (int i = 1; i <= mesh.TotalPoints; i++)
            {
                mesh.points[i].x = mesh.points[i].x + U[DOF*i-2];
                mesh.points[i].y = mesh.points[i].y + U[DOF*i-1];
                mesh.points[i].z = mesh.points[i].z + U[DOF*i];
                deformation[i] = sqrt(U[DOF*i-2]*U[DOF*i-2] + U[DOF*i-1]*U[DOF*i-1] + U[DOF*i]*U[DOF*i]);
            }
        }break;
        case 2:
        {
            int load_type = 1;
            vector<double> U_supp(DOF*mesh.SurfaceNds+1);
            vector<double> dof(DOF*mesh.SurfaceNds+1);
            vector<vector<double> > D_matrix = computeD2D(load_type, materials.MOE, materials.nu, materials.thickness);
            
            for (int g = 1; g <= TotQuadPnts; g++)
            {
                int e = kernels[g].elem_num;
                vector<vector<double> > B = computeB2D(mesh.SurfaceNds, kernels[g].dNdx, kernels[g].dNdy);
                int kk = 1;
                for (int j = 1; j <= mesh.SurfaceNds; j++)
                {
                    for (int iter = DOF-1; iter >= 0; iter--)
                    {
                        dof[kk] =  DOF*mesh.surfaces[e][j]-iter;
                        kk = kk + 1;
                    }
                }		
                for (int i = 1; i <= DOF*mesh.SurfaceNds; i++)
                {U_supp[i] = U[dof[i]];}
            
                double strain_xx = 0;
                double strain_yy = 0;
                double gamma_xy = 0;
                for (int s = 1; s <= DOF*mesh.SurfaceNds; s++)
                {
                    strain_xx = strain_xx  + B[s][1]*U_supp[s];
                    strain_yy = strain_yy  + B[s][2]*U_supp[s];
                    gamma_xy = gamma_xy  + B[s][3]*U_supp[s];
                }
                double stress_xx = D_matrix[1][1]*strain_xx + D_matrix[2][1]*strain_yy + D_matrix[3][1]*gamma_xy;
                double stress_yy = D_matrix[1][2]*strain_xx + D_matrix[2][2]*strain_yy + D_matrix[3][2]*gamma_xy;
                //stress_xy = D_matrix[1][3]*strain_xx + D_matrix[2][3]*strain_yy + D_matrix[3][3]*gamma_xy;
                gauss_vms[g] = sqrt(pow((stress_xx-stress_yy),2)+ stress_xx*stress_xx + stress_yy*stress_yy + 6*gamma_xy*gamma_xy)/sqrt(2.0);
            }
        
            cout << "Computing final nodal deformation..." << endl;
            for (int i = 1; i <= mesh.TotalPoints; i++)
            {
                mesh.points[i].x = mesh.points[i].x + U[DOF*i-1];
                mesh.points[i].y = mesh.points[i].y + U[DOF*i];
                deformation[i] = sqrt(U[DOF*i-1]*U[DOF*i-1] + U[DOF*i]*U[DOF*i]);
            }
        }break;
    }
	//interpolate stresses from gauss points to nodes
	for (int i = 1; i <= mesh.TotalPoints; i++)
	{
		double sum_w = 0;
		double sum_vms = 0;
		for (int j = 1; j <= ndkernels[i].TotalNbrs; j++)
		{
			int e = ndkernels[i].nbrs[j];
			for (int g = 1; g <= elems[e].TotalNbrs; g++)
			{
				int q = elems[e].nbrs[g];
				double dx = mesh.points[i].x - kernels[q].pnt.x;
				double dy = mesh.points[i].y - kernels[q].pnt.y;
				double dz = mesh.points[i].z - kernels[q].pnt.z;
				double d = sqrt(dx*dx+dy*dy+dz*dz);
				double w = 1/d;
				sum_w = sum_w + w;
				sum_vms = sum_vms + w*gauss_vms[q];
			}
		}
		vms[i] = sum_vms/sum_w;
		FOS[i] = vms[i]/materials.YTS;
	}

	cout << "Printing results to file..." << endl;
	mesh.printMeshVTK(deformation, "U", 1);
	mesh.printMeshVTK(vms,"vms",1);
	mesh.printMeshVTK(FOS,"fos",1);
	
	mesh.printTXT(U, "U", 1);
    mesh.printTXT(vms, "vms", 1);
    mesh.printTXT(FOS, "fos", 1);

    cout << "SOLVER COMPLETE" << endl;   
}
