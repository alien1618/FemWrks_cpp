#include "../femwrks.h"

void GLOBAL_MATRICES::assembleDiffusionMatrix(bool SUPG, int TotalMeshPoints, vector<vector<int>> elements, int elem_nds, vector<KERNEL> interps, int gauss_TotalPoints, double D, vector<double> Vx1, vector<double> Vy1, vector<double> Vz1, double stab_prmtr)
{
    cout << "SUPG = " << SUPG << endl;
    K.resize(TotalMeshPoints+1, vector<double> (TotalMeshPoints+1));
    cout << "Assembling diffusion matrix..." << endl;
    for (int g = 1; g <= gauss_TotalPoints; g++)
    {
        int e = interps[g].elem_num;
        int cmsize = elem_nds+1;
        double Vx = 0;
        double Vy = 0;
        double Vz = 0;
        for (int mm = 1; mm <= elem_nds; mm++)
        {
            int cm = elements[e][mm];
            Vx = Vx + interps[g].N[mm]*Vx1[cm];
            Vy = Vy + interps[g].N[mm]*Vy1[cm];
            Vz = Vz + interps[g].N[mm]*Vz1[cm];
        }

        vector<vector<double> > K_gauss_diffusion(cmsize,vector<double> (cmsize));
        vector<vector<double> > K_gauss_advection(cmsize,vector<double> (cmsize));

        for (int j = 1 ; j <= elem_nds; j++)
        {
            for (int i = 1; i <= elem_nds; i++)
            {
                K_gauss_diffusion[i][j] = D*(interps[g].dNdNx[i][j] + interps[g].dNdNy[i][j] + interps[g].dNdNz[i][j]);
                K_gauss_advection[i][j] = interps[g].NdNx[i][j]*Vx+interps[g].NdNy[i][j]*Vy+interps[g].NdNz[i][j]*Vz;
            }
        }

		// cout << "ASSEMBLE THE GAUSS STIFFNESS AND FORCE MATRICES DIRECTLY INTO THE GLOBAL MATRICES" << endl;
        for (int n = 1; n <= elem_nds; n++)
        {
            for (int m = 1; m <= elem_nds; m++)
            {
                int cm = elements[e][m];
                int cn = elements[e][n];
                K[cm][cn] = K[cm][cn] + interps[g].quad_weight * (K_gauss_diffusion[m][n]+K_gauss_advection[m][n]);
            }
        }
        if (SUPG == true)
        {
            //cout << "COMPUTE GAUSS MATRICES FOR STABILIZATION SCHEMES" << endl;
            double norm_V = pow((Vx*Vx + Vy*Vy + Vz*Vz),0.5);
            double tau_SUPG = 0;
            if (norm_V == 0)
            {tau_SUPG = (0.5*stab_prmtr);}
            else
            {tau_SUPG = (stab_prmtr)/norm_V;}
            
            vector<double> dN_V(elem_nds+1);
			vector<vector<double> > K_gauss_SUPG(elem_nds+1,vector<double> (elem_nds+1));
    		for (int n = 1; n <= elem_nds; n++)
    		{	
        		for (int m = 1; m <= elem_nds; m++)
        		{	
            		K_gauss_SUPG[m][n] = 0;
        		}
    		}
			for (int n = 1;  n <= elem_nds; n++)
			{
                int cj = elements[e][n];
                dN_V[n] = Vx1[cj]*interps[g].dNdx[n] + Vy1[cj]*interps[g].dNdy[n] + Vz1[cj]*interps[g].dNdz[n];
            }
			K_gauss_SUPG = multiply_transp(dN_V, dN_V, elem_nds);
			//M_gauss_SUPG = multiply_transp(dN_V, interps[g].N, elem_nds);
			//for (int n = 1; n <= elem_nds; n++)
			//{F_gauss_SUPG[n] = dN_V[n]*Q;}
			
            //cout << "ASSEMBLE THE GAUSS MATRICES DIRECTLY INTO THE GLOBAL MATRICES" << endl;
            for (int n = 1; n <= elem_nds; n++)
            {
                for (int m = 1; m <= elem_nds; m++)
                {
                    int cm = elements[e][m];
                    int cn = elements[e][n];
                    K[cm][cn] = K[cm][cn] + interps[g].quad_weight *(tau_SUPG*K_gauss_SUPG[m][n]);
                }
            }
        }
    }
}
void GLOBAL_MATRICES::assembleMassMatrix(bool SUPG, int TotalMeshPoints, vector<vector<int>> elements, int elem_nds, vector<KERNEL> interps, int gauss_TotalPoints, vector<double> Vx1, vector<double> Vy1, vector<double> Vz1, double stab_prmtr)
{
    cout << "SUPG = " << SUPG << endl;
    M.resize(TotalMeshPoints+1, vector<double> (TotalMeshPoints+1));
    cout << "Assembling mass matrix..." << endl;
    for (int g = 1; g <= gauss_TotalPoints; g++)
    {
        int e = interps[g].elem_num;
        int cmsize = elem_nds+1;
        double Vx = 0;
        double Vy = 0;
        double Vz = 0;
        for (int mm = 1; mm <= elem_nds; mm++)
        {
            int cm = elements[e][mm];
            Vx = Vx + interps[g].N[mm]*Vx1[cm];
            Vy = Vy + interps[g].N[mm]*Vy1[cm];
            Vz = Vz + interps[g].N[mm]*Vz1[cm];
        }
        vector<vector<double> > M_gauss(cmsize,vector<double> (cmsize));
        for (int j = 1 ; j <= elem_nds; j++)
        {
            for (int i = 1; i <= elem_nds; i++)
            {
                M_gauss[i][j] = interps[g].NN[i][j];
            }
        }

        // cout << "COMPUTE GAUSS MASS MATRICES" << endl;
        for (int j = 1 ; j <= elem_nds; j++)
        {
            for (int i = 1; i <= elem_nds; i++)
            {M_gauss[i][j] = interps[g].NN[i][j];}
        }
        // cout << "ASSEMBLE THE GAUSS MASS MATRIX DIRECTLY INTO THE GLOBAL MASS MATRIX" << endl;
        for (int n = 1; n <= elem_nds; n++)
        {
            for (int m = 1; m <= elem_nds; m++)
            {
                int cm = elements[e][m];
                int cn = elements[e][n];
                M[cm][cn] = M[cm][cn] + interps[g].quad_weight * (M_gauss[m][n]);
            }
        }

        if (SUPG == true)
        {
            //cout << "COMPUTE GAUSS MATRICES FOR STABILIZATION SCHEMES" << endl;
            double norm_V = pow((Vx*Vx + Vy*Vy + Vz*Vz),0.5);
            double tau_SUPG = 0;
            if (norm_V == 0)
            {tau_SUPG = (0.5*stab_prmtr);}
            else
            {tau_SUPG = (stab_prmtr)/norm_V;}
            
            vector<double> dN_V(elem_nds+1);
			vector<vector<double> > M_gauss_SUPG(elem_nds+1,vector<double> (elem_nds+1));
    		for (int n = 1; n <= elem_nds; n++)
    		{	
        		for (int m = 1; m <= elem_nds; m++)
        		{	
            		M_gauss_SUPG[m][n] = 0;
        		}
    		}
			for (int n = 1;  n <= elem_nds; n++)
			{
                int cj = elements[e][n];
                dN_V[n] = Vx1[cj]*interps[g].dNdx[n] + Vy1[cj]*interps[g].dNdy[n] + Vz1[cj]*interps[g].dNdz[n];
            }
			//M_gauss_SUPG = multiply_transp(dN_V, interps[g].N, elem_nds);
			//for (int n = 1; n <= elem_nds; n++)
			//{F_gauss_SUPG[n] = dN_V[n]*Q;}
			
            //cout << "ASSEMBLE THE GAUSS MATRICES DIRECTLY INTO THE GLOBAL MATRICES" << endl;
            for (int n = 1; n <= elem_nds; n++)
            {
                for (int m = 1; m <= elem_nds; m++)
                {
                    int cm = elements[e][m];
                    int cn = elements[e][n];
                    M[cm][cn] = M[cm][cn] + interps[g].quad_weight*tau_SUPG*M_gauss_SUPG[m][n];
                }
            }
        }
    }
}
void GLOBAL_MATRICES::assembleStiffnessMatrix(int dim, int load_type, double E, double nu, double t, int TotalNds, vector<KERNEL> kernels, int totquadpnts,  vector<vector<int> > elements, int elemnds)
{
	int DOF = dim;
    double msize = DOF*TotalNds+1;
    K.resize(msize, vector<double> (msize));
    F.resize(msize);
    for (int s = 1; s <= totquadpnts; s++)
    {
        int e = kernels[s].elem_num;
        switch(dim)
        {
            case 2:
            {
                vector<vector<double> > B(DOF*elemnds+1, vector<double>(4));
                vector<double> dof(DOF*elemnds+1);
                vector<vector<double> > B_transp(4, vector<double>(DOF*elemnds+1));
                vector<vector<double> > temp_mat10(4, vector<double>(DOF*elemnds+1));
                vector<vector<double> > D_matrix(4, vector<double>(4));
                vector<vector<double> > K_gauss_diffusion(DOF*elemnds+1,vector<double> (DOF*elemnds+1));

                //cout << "compute gauss matrices" << endl;
                D_matrix = computeD2D(load_type,E,nu,t);
                B = computeB2D(elemnds, kernels[s].dNdx, kernels[s].dNdy);
                B_transp = transpose(B, DOF*elemnds, 3);
                temp_mat10 = multiply(B_transp, 3, DOF*elemnds, D_matrix, 3, 3);
                K_gauss_diffusion = multiply(temp_mat10, 3, DOF*elemnds, B, DOF*elemnds, 3);
                //cout << "ASSEMBLE THE GAUSS MATRICES DIRECTLY INTO THE GLOBAL MATRICES" << endl;
                
                int kk = 1;
                for (int j = 1; j <= elemnds; j++)
                {
                    for (int iter = DOF-1; iter >= 0; iter--)
                    {
                        dof[kk] =  DOF*elements[e][j]-iter;
                        kk = kk + 1;
                    }
                }
                for (int l = 1; l <= DOF*elemnds; l++)
                {
                    for (int j = 1; j <= DOF*elemnds; j++)
                    {K[dof[j]][dof[l]] = K[dof[j]][dof[l]] + K_gauss_diffusion[l][j]*kernels[s].quad_weight;}
                }
                
            }break;
            case 3:
            {
                vector<vector<double> > B(DOF*elemnds+1, vector<double>(7));
                vector<double> dof(DOF*elemnds+1);
                vector<vector<double> > B_transp(7, vector<double>(DOF*elemnds+1));
                vector<vector<double> > temp_mat10(7, vector<double>(DOF*elemnds+1));
                vector<vector<double> > D_matrix(7, vector<double>(7));
                vector<vector<double> > K_gauss_diffusion(DOF*elemnds+1,vector<double> (DOF*elemnds+1));

                //cout << "compute gauss matrices" << endl;
                D_matrix = computeD3D(E,nu);
                B = computeB3D(elemnds, kernels[s].dNdx, kernels[s].dNdy, kernels[s].dNdz);
                B_transp = transpose(B, DOF*elemnds, 6);
                temp_mat10 = multiply(B_transp, 6, DOF*elemnds, D_matrix, 6, 6);
                K_gauss_diffusion = multiply(temp_mat10, 6, DOF*elemnds, B, DOF*elemnds, 6);
                //cout << "ASSEMBLE THE GAUSS MATRICES DIRECTLY INTO THE GLOBAL MATRICES" << endl;
                int kk = 1;
                for (int j = 1; j <= elemnds; j++)
                {
                    for (int iter = DOF-1; iter >= 0; iter--)
                    {
                        dof[kk] =  DOF*elements[e][j]-iter;
                        kk = kk + 1;
                    }
                }
                for (int l = 1; l <= DOF*elemnds; l++)
                {
                    for (int j = 1; j <= DOF*elemnds; j++)
                    {K[dof[j]][dof[l]] = K[dof[j]][dof[l]] + K_gauss_diffusion[l][j]*kernels[s].quad_weight;}
                }
                    
            }break;
        }
    }
}
void GLOBAL_MATRICES::assembleForceMatrix(bool SUPG, int TotalMeshPoints, vector<vector<int>> elements, int elem_nds, vector<KERNEL> interps, int gauss_TotalPoints, vector<double> Vx1, vector<double> Vy1, vector<double> Vz1, vector<double> Q1, double stab_prmtr)
{
    cout << "SUPG = " << SUPG << endl;
    F.resize(TotalMeshPoints+1);
    cout << "Assembling global matrices..." << endl;
    for (int g = 1; g <= gauss_TotalPoints; g++)
    {
        int e = interps[g].elem_num;
        int cmsize = elem_nds+1;
        double Vx = 0;
        double Vy = 0;
        double Vz = 0;
        double Q = 0;
        for (int mm = 1; mm <= elem_nds; mm++)
        {
            int cm = elements[e][mm];
            Vx = Vx + interps[g].N[mm]*Vx1[cm];
            Vy = Vy + interps[g].N[mm]*Vy1[cm];
            Vz = Vz + interps[g].N[mm]*Vz1[cm];
            Q = Q + interps[g].N[mm]*Q1[cm];
        }

        vector<double>  F_gauss(cmsize);

        for (int j = 1 ; j <= elem_nds; j++)
        {
             F_gauss[j] =  interps[g].N[j]*Q;
        }

		// cout << "ASSEMBLE THE GAUSS STIFFNESS AND FORCE MATRICES DIRECTLY INTO THE GLOBAL MATRICES" << endl;
        for (int n = 1; n <= elem_nds; n++)
        {
            int co = elements[e][n];
            F[co] = F[co] + interps[g].quad_weight * (F_gauss[n]);
        }
        if (SUPG == true)
        {
            //cout << "COMPUTE GAUSS MATRICES FOR STABILIZATION SCHEMES" << endl;
            /*
            double norm_V = pow((Vx*Vx + Vy*Vy + Vz*Vz),0.5);
            double tau_SUPG = 0;
            if (norm_V == 0)
            {tau_SUPG = (0.5*stab_prmtr);}
            else
            {tau_SUPG = (stab_prmtr)/norm_V;}
            */
            vector<double> dN_V(elem_nds+1);
			vector<double>  F_gauss_SUPG(elem_nds+1);
    		for (int n = 1; n <= elem_nds; n++)
    		{	
        		F_gauss_SUPG[n] = 0;
    		}
			for (int n = 1;  n <= elem_nds; n++)
			{
                int cj = elements[e][n];
                dN_V[n] = Vx1[cj]*interps[g].dNdx[n] + Vy1[cj]*interps[g].dNdy[n] + Vz1[cj]*interps[g].dNdz[n];
            }
			//for (int n = 1; n <= elem_nds; n++)
			//{F_gauss_SUPG[n] = dN_V[n]*Q;}
			
            //cout << "ASSEMBLE THE GAUSS MATRICES DIRECTLY INTO THE GLOBAL MATRICES" << endl;
            for (int n = 1; n <= elem_nds; n++)
            {
                int co = elements[e][n];
                F[co] = F[co] + interps[g].quad_weight *F_gauss_SUPG[n];
            }
        }
    }
}
void GLOBAL_MATRICES::applyDBCDirectTransient(int TotalPoints, vector<int> bc_pnts, vector<double> bc_val, int bc_totpnts)
{
    if (bc_totpnts > 0)
    {
        for (int k = 1; k <= bc_totpnts; k++)
        {
            for (int i = 1; i <= TotalPoints; i++)
            {F[i] = F[i]- (K[bc_pnts[k]][i]*bc_val[k]);}
        }
        for (int h = 1; h <=bc_totpnts; h++)
        {F[bc_pnts[h]] = bc_val[h];}
        for (int h = 1; h <= bc_totpnts; h++)
        {
            for (int k = 1; k <= TotalPoints; k++)
            {
                K[bc_pnts[h]][k] = 0;
                K[k][bc_pnts[h]] = 0;
                M[bc_pnts[h]][k] = 0;
                M[k][bc_pnts[h]] = 0;
            }
            K[bc_pnts[h]][bc_pnts[h]] = 1;
            M[bc_pnts[h]][bc_pnts[h]] = 1;
        }
    }
}
void GLOBAL_MATRICES::applyDBCDirectSteadyState(int DOF, int TotalPoints, vector<int> bc_pnts, vector<double> bc_val, int bc_totpnts)
{
    if (bc_totpnts > 0)
    {
        for (int k = 1; k <= bc_totpnts; k++)
        {
            for (int m = 1; m <= DOF; m++)
            {
                int iter = DOF-m;
                for (int i = 1; i <= DOF*TotalPoints; i++)
                {
                    //cout << "values " << values[DOF*k-iter] << endl;
                    //F[i] = F[i]- (K[DOF*bc_pnts[k]-iter][i]*values[DOF*k-iter]);
                    F[i] = F[i]-(K[DOF*bc_pnts[k]-iter][i]*bc_val[k]);
                }
            }
        }
        for (int h = 1; h <=bc_totpnts; h++)
        {
            for (int m = 1; m <= DOF; m++)
            {
                int iter = DOF-m;
               // cout << "DOF*bc_pnts[h]-iter " << DOF*bc_pnts[h]-iter << " val " << bc_val[h] << endl;
                //F[DOF*bc_pnts[h]-iter] = values[DOF*h-iter];
                F[DOF*bc_pnts[h]-iter] = bc_val[h];
            }
        }
        for (int h = 1; h <= bc_totpnts; h++)
        {
            for (int m = 1; m <= DOF; m++)
            {
                int iter = DOF-m;
                for (int k = 1; k <= DOF*TotalPoints; k++)
                {
                    K[DOF*bc_pnts[h]-iter][k] = 0;
                    K[k][DOF*bc_pnts[h]-iter] = 0;
                }
                K[DOF*bc_pnts[h]-iter][DOF*bc_pnts[h]-iter] = 1;
            }
        }
    }
}
void GLOBAL_MATRICES::applyDBCPenaltyTransient(int TotalPoints, vector<int> bc_pnts, vector<double> bc_val, int bc_totpnts)
{
    if (bc_totpnts > 0)
    {
        vector<double> diagonal(TotalPoints+1);
        for (int i = 1; i <= TotalPoints; i++)
        {diagonal[i] = K[i][i];}
        double maxK = maximum <double> (diagonal, TotalPoints);
        double alpha = maxK*10E8;

        for (int h = 1; h <= bc_totpnts; h++)
        {F[bc_pnts[h]] = alpha*K[bc_pnts[h]][bc_pnts[h]]*bc_val[h];}
        for (int h = 1; h <= bc_totpnts; h++)
        {
            double tmp1 = K[bc_pnts[h]][bc_pnts[h]];
            double tmp2 = M[bc_pnts[h]][bc_pnts[h]];
            K[bc_pnts[h]][bc_pnts[h]] = alpha*tmp1;
            M[bc_pnts[h]][bc_pnts[h]] = alpha*tmp2;
        }
    }
}
GLOBAL_MATRICES::~GLOBAL_MATRICES ()
{
    //destructor
}
vector<double> advanceInTime(int transient, vector<double> U_now, GLOBAL_MATRICES global_matrices, int TotalPoints, double dt, int CGiter,  double tol)
{
    cout << "Solving the system of equations..." << endl;
    vector<double> U_future(TotalPoints+1);
    switch(transient)
    {
        case 0:
        {U_future = solvePCG(global_matrices.K,  global_matrices.F, TotalPoints, CGiter, tol);}break;
        case 1:
        {U_future = solveTransientImplicitPCG(global_matrices.K, global_matrices.M, global_matrices.F, U_now, TotalPoints, dt, CGiter, tol);}break;
        case 2:
        {U_future = solveTransientExplicitPCG(global_matrices.K, global_matrices.M, global_matrices.F, U_now, TotalPoints, dt, CGiter, tol);}break;
    }
    return U_future;
}
