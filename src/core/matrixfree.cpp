#include "../femwrks.h"

vector<NBR> computeNodeKernels(int totnds, vector<vector<int>> elements, int TotElems, int elemnds)
{
    vector<NBR> ndkernels(totnds+1);
    for (int i = 1; i<= totnds; i++)
    {
        ndkernels[i].TotalNbrs = 0;
        ndkernels[i].nbrs.resize(1);
    }
    for (int j = 1; j <= TotElems; j++)
    {
        for (int k = 1; k <= elemnds; k++)
        {
            int ndnum = elements[j][k];
            ndkernels[ndnum].TotalNbrs++;
            ndkernels[ndnum].nbrs.push_back(1);
            ndkernels[ndnum].nbrs[ndkernels[ndnum].TotalNbrs] = j;
        }
    }
    return ndkernels;
}

tuple <vector<int>, int, vector<int>, int> collectGaussPointsWithBC(vector<int> bc_ndnum, int bc_totnds, vector<vector<int>> elements, int elemnds, vector<KERNEL> kernels, int TotalGaussPoints)
{
    int tot_g_bc = 0;
    vector<int> g_bc(1);
    vector<int> indx(TotalGaussPoints+1);
    for (int g = 1; g <= TotalGaussPoints; g++)
    {indx[g] = 0;}
    for (int g = 1; g <= TotalGaussPoints; g++)
    {
        int e = kernels[g].elem_num;
        for (int k = 1; k <= bc_totnds; k++)
        {
        //  for (int i = 1; i <= mesh.TotalPoints; i++)
        //  {
            
            for (int n = 1; n <= elemnds; n++)
            {
                for (int m = 1; m <= elemnds; m++)
                {
                    int cm = elements[e][m];
                    //int cn = elements[e].points[n];
                    if (cm == bc_ndnum[k] && indx[g] == 0)
                    //if (cm == bc_ndnum[k] && cn == i && indx[g] == 0)
                    {
                        tot_g_bc++;
                        g_bc.push_back(1);
                        g_bc[tot_g_bc] = g;
                        indx[g] = 1;
                    }
                 }
             }
        // }
        }
    }
    int tot_g_bc2 = 0;
    vector<int> g_bc2(1);
    vector<int> indx2(TotalGaussPoints+1);
    for (int g = 1; g <= TotalGaussPoints; g++)
    {indx2[g] = 0;}
    for (int g = 1; g <= TotalGaussPoints; g++)
    {
        int e = kernels[g].elem_num;
        for(int j = 1; j <= elemnds; j++)
        {
           for(int i = 1; i <= elemnds; i++)
           {
                int ci = elements[e][i];
                int check = 0;
                for (int h = 1; h <= bc_totnds; h++)
                {
                	if (ci == bc_ndnum[h])
                	{check++;}
                }
                if (check == 0 && indx2[g] == 0)
                {
                    tot_g_bc2++;
                    g_bc2.push_back(1);
                    g_bc2[tot_g_bc2] = g;
                    indx2[g] = 1;
                }
            }
        }
    }
    return make_tuple(g_bc, tot_g_bc, g_bc2, tot_g_bc2);
}
vector<double> solveElementByElement(vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> Q, double D, int totnds, vector<KERNEL> kernels, int TotalGaussPoints, vector<vector<int>> elements, int elemnds, vector<NBR> ndkernels, double dt, int transient, bool SUPG, double stab_prmtr, int CGiter, double conv_tol, vector<int> bc_ndnum, vector<double> bc_ndval, int bc_totnds, vector<int> g_bc, int tot_g_bc, vector<int> g_bc2, int tot_g_bc2)
{
	//this function solves the steady state or transient solution of
	//the general transport equation using the finite element method
	//without generating global matrices.
	//it solves it element-by-element iteratively using the pre-conditioned
	//conjugate gradient method. Function assumes there is only one
	//degree of freedom per node. Used for heat/fluid flow simulations.
	
    vector<double> F(totnds+1);
    for(int j = 1; j <= totnds; j++)
    {F[j] = 0;}
    vector<vector<double> > MK(elemnds+1,vector<double> (elemnds+1));
    for (int g = 1; g <= TotalGaussPoints; g++)
    {
        //----------------------------------------------------------
        //cout << "interpolate material properties at gauss point" << endl;
        //----------------------------------------------------------
        double Vx_g = 0;
        double Vy_g = 0;
        double Vz_g = 0;
        double Q_g = 0;
        int e = kernels[g].elem_num;
        for (int mm = 1; mm <= elemnds; mm++)
        {
            int cm = elements[e][mm];
            Vx_g = Vx_g + kernels[g].N[mm]*Vx[cm];
            Vy_g = Vy_g + kernels[g].N[mm]*Vy[cm];
            Vz_g = Vz_g + kernels[g].N[mm]*Vz[cm];
            Q_g = Q_g + kernels[g].N[mm]*Q[cm];
        }

        //----------------------------------------------------------
        // cout << "COMPUTE GAUSS STIFFNESS AND FORCE MATRICES" << endl;
        //----------------------------------------------------------
        if (transient == 0)
        {
            //CALCULATES IMPLICIT FORM: 0 = (K*U + Q)
            //First: calculate K
            for (int j = 1 ; j <= elemnds; j++)
            {
                for (int i = 1; i <= elemnds; i++)
                {
                    double advection_term = kernels[g].NdNx[i][j]*Vx_g + kernels[g].NdNy[i][j]*Vy_g+ kernels[g].NdNz[i][j]*Vz_g;
                    double diffusion_term = D*(kernels[g].dNdNx[i][j]+kernels[g].dNdNy[i][j]+kernels[g].dNdNz[i][j]);
                    MK[i][j] = kernels[g].quad_weight*(diffusion_term+advection_term);
                }
            }
            kernels[g].MK = MK;

            //Second: calculate (Q)
            for (int j = 1 ; j <= elemnds; j++)
            {
                int cj = elements[e][j];
                F[cj] = F[cj] + kernels[g].quad_weight*(dt*Q_g);
            }
        }
        if (transient == 1)
        {
			//CALCULATES IMPLICIT FORM: U_new(M + K*dt) = ((M*U_now) + (Q*dt))
			//First: calculate (M + K*dt)
			for (int j = 1 ; j <= elemnds; j++)
			{
				for (int i = 1; i <= elemnds; i++)
				{
					int cj = elements[e][i];
					double advection_term = kernels[g].NdNx[i][j]*Vx[cj] + kernels[g].NdNy[i][j]*Vy[cj] + kernels[g].NdNz[i][j]*Vz[cj];
					double diffusion_term = D*(kernels[g].dNdNx[i][j]+ kernels[g].dNdNy[i][j]+ kernels[g].dNdNz[i][j]);
					double capacitance_term = kernels[g].NN[i][j];
					MK[i][j] = kernels[g].quad_weight*(capacitance_term+dt*(diffusion_term+advection_term));
				}
			}
			if (SUPG == true)
			{
				double norm_V = pow((Vx_g*Vx_g + Vy_g*Vy_g + Vz_g*Vz_g),0.5)+1e-5;
				double tau_SUPG = stab_prmtr/norm_V;
				vector<double> dN_V(elemnds+1);
				vector<vector<double> > K_gauss_SUPG(elemnds+1,vector<double> (elemnds+1));
				//vector<vector<double> > M_gauss_SUPG(elemnds+1,vector<double> (elemnds+1));

				for (int n = 1;  n <= elemnds; n++)
				{
					int cj = elements[e][n];
					dN_V[n] = Vx[cj]*kernels[g].dNdx[n] + Vy[cj]*kernels[g].dNdy[n] + Vz[cj]*kernels[g].dNdz[n];
				}
				K_gauss_SUPG = multiply_transp(dN_V, dN_V, elemnds);
				for (int j = 1 ; j <= elemnds; j++)
				{
					for (int i = 1; i <= elemnds; i++)
					{
						MK[i][j] = MK[i][j] + kernels[g].quad_weight*dt*tau_SUPG*K_gauss_SUPG[i][j];
					}
				}
			}
			kernels[g].MK = MK;
			//for (int n = 1; n <= c_totnds; n++)
			//{F_gauss_SUPG[n] = dN_V[n]*Q;}
			//Second: calculate ((M * U_now) + (Q*dt))
			for (int j = 1 ; j <= elemnds; j++)
			{
				int cj = elements[e][j];
				double row_sum = 0;
				for (int i = 1; i <= elemnds; i++)
				{
					int ci = elements[e][i];
					row_sum = row_sum + (kernels[g].NN[i][j])*U[ci];
				}
				F[cj] = F[cj] + kernels[g].quad_weight*(row_sum + dt*Q_g);
			}
        }
    }

    U = solvePCG(1, F, ndkernels, totnds, elements, elemnds, kernels, bc_ndnum, bc_ndval, bc_totnds, g_bc, tot_g_bc, g_bc2, tot_g_bc2, CGiter, conv_tol);
    return U;
}

vector<double> solveElementByElement(int dim, int DOF, double E, double nu, double t, int totnds, vector<KERNEL> kernels, int TotalGaussPoints, vector<vector<int>> elements, int elemnds, vector<NBR> ndkernels, int CGiter, double conv_tol, vector<double> F, vector<int> bc_ndnum, vector<double> bc_ndval, int bc_totnds, vector<int> g_bc, int tot_g_bc, vector<int> g_bc2, int tot_g_bc2)
{
	//this function solves the elasticity equation using the finite 
	//element method without generating global matrices.
	//it solves it element-by-element iteratively using the pre-conditioned
	//conjugate gradient method. Function assumes either 2 DOF for 2D or 3 DOF
	// for 3D per node. Used for 2D and 3D linear stress analysis.
	
    vector<vector<double> > MK(DOF*elemnds+1,vector<double> (DOF*elemnds+1));
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

			for (int g = 1; g <= TotalGaussPoints; g++)
			{
				//cout << "compute gauss matrices" << endl;
				D_matrix = computeD2D(1,E,nu,t);
				B = computeB2D(elemnds, kernels[g].dNdx, kernels[g].dNdy);
				B_transp = transpose(B, DOF*elemnds, 3);
				temp_mat10 = multiply(B_transp, 3, DOF*elemnds, D_matrix, 3, 3);
				K_gauss_diffusion = multiply(temp_mat10, 3, DOF*elemnds, B, DOF*elemnds, 3);
				   
				//First: calculate K
				for (int j = 1 ; j <= DOF*elemnds; j++)
				{
					for (int i = 1; i <= DOF*elemnds; i++)
					{
						MK[i][j] = kernels[g].quad_weight*K_gauss_diffusion[i][j];
					}
				}
				kernels[g].MK = MK;
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

			for (int g = 1; g <= TotalGaussPoints; g++)
			{

				//cout << "compute gauss matrices" << endl;
				D_matrix = computeD3D(E,nu);
				B = computeB3D(elemnds, kernels[g].dNdx, kernels[g].dNdy, kernels[g].dNdz);
				B_transp = transpose(B, DOF*elemnds, 6);
				temp_mat10 = multiply(B_transp, 6, DOF*elemnds, D_matrix, 6, 6);
				K_gauss_diffusion = multiply(temp_mat10, 6, DOF*elemnds, B, DOF*elemnds, 6);
			
				//First: calculate K
				for (int j = 1 ; j <= DOF*elemnds; j++)
				{
					for (int i = 1; i <= DOF*elemnds; i++)
					{
						MK[i][j] = kernels[g].quad_weight*K_gauss_diffusion[i][j];
					}
				}
				kernels[g].MK = MK;
			}
		}break;
	}
    vector<double> U = solvePCG(DOF, F, ndkernels, totnds, elements, elemnds, kernels, bc_ndnum, bc_ndval, bc_totnds, g_bc, tot_g_bc, g_bc2, tot_g_bc2, CGiter, conv_tol);
    return U;
}

vector<double> solvePCG(int DOF, vector<double> F, vector<NBR> ndkernels, int totnds, vector<vector<int>> elements, int elemnds, vector<KERNEL> kernels, vector<int> bc_ndnum, vector<double> bc_ndval, int bc_totnds, vector<int> g_bc, int tot_g_bc, vector<int> g_bc2, int tot_g_bc2, int CGiter, double conv_tol)
{
    /*
    // applying boundary conditions to the force vector
    for (int k = 1; k <= bc_totnds; k++)
    {
        for (int i = 1; i <= totnds; i++)
        {
        int ndnum = bc_ndnum[k];
        for (int s = 1; s <= ndkernels[ndnum].TotalElems; s++)
        {
            int g = ndkernels[ndnum].elems[s]; // careful...assumes a single gauss point per element
            double sum = 0;
            for (int n = 1; n <= elemnds; n++)
            {
            for (int m = 1; m <= elemnds; m++)
            {
                int cm = elements[e].points[m];
                int cn = elements[e].points[n];
                if (cm == ndnum && cn == i)
                {sum = sum + kernels[g].MK[m][n]*bc_ndval[k];}
            }
            }
            F[i] = F[i]-sum;
        }
        }
    }
    */

    // applying boundary conditions to the force vector	
    for (int s = 1; s <= tot_g_bc; s++)
    {
        int g = g_bc[s];
        int e = kernels[g].elem_num;
                
        for (int k = 1; k <= bc_totnds; k++)
        {
			for(int q = 1; q <= DOF; q++)
			{
				int iter = DOF-q;
				for (int i = 1; i <= totnds; i++)
				{
					double sum = 0;
					for (int n = 1; n <= elemnds; n++)
					{
						for (int m = 1; m <= elemnds; m++)
						{
							int cm = elements[e][m];
							int cn = elements[e][n];
							if (cm == bc_ndnum[k] && cn == i)
							{sum = sum + kernels[g].MK[DOF*m-iter][DOF*n-iter]*bc_ndval[k];}
						}
					}
					F[DOF*i-iter] = F[DOF*i-iter]-sum;
				}
			}
        }
    }

    for (int h = 1; h <= bc_totnds; h++)
    {
        for (int m = 1; m <= DOF; m++)
        {
            int iter = DOF-m;
            F[DOF*bc_ndnum[h]-iter] = bc_ndval[h];
        }
    }

    //Conjugate Gradient Solver(Serial Code with no OpenMP)
    vector<double> U_new(DOF*totnds+1);
    vector<double> U_new_iter(DOF*totnds+1);
    vector<double> A_x(DOF*totnds+1);
    vector<double> A_p(DOF*totnds+1);
    vector<double> r(DOF*totnds+1);
    vector<double> r_new(DOF*totnds+1);
    vector<double> p(DOF*totnds+1);

    U_new = F; // initial guess of U_new
    for (int s = 1; s <= tot_g_bc2; s++)
    {
        int g = g_bc2[s];
        int e = kernels[g].elem_num;
		
        vector<int> dof(DOF*elemnds+1);
        int kk = 1;
        for (int f = 1; f <= elemnds; f++)
        {
            for (int iter = DOF-1; iter >= 0; iter--)
            {
                dof[kk] =  DOF*elements[e][f]-iter;
                kk = kk + 1;
            }
        }
        for(int j = 1; j <= DOF*elemnds; j++)
        {    
            double product_gauss = 0; 
            for(int i = 1; i <= DOF*elemnds; i++)
            {
                if (DOF == 1)
                {
                    int ci = elements[e][i];
                    int check = 0;
                
                    for (int h = 1; h <= bc_totnds; h++)
                    {
                        if (ci == bc_ndnum[h])
                        {check++;}
                    }
                    if (check == 0)
                    {product_gauss = product_gauss + kernels[g].MK[i][j]*U_new[ci];}
                }
                else
                {product_gauss = product_gauss + kernels[g].MK[i][j]*U_new[dof[i]];}
            }
            A_x[dof[j]] = A_x[dof[j]] + product_gauss;
        }
    }
    for (int k = 1; k <= bc_totnds; k++)
    {
		for (int m = 1; m <= DOF; m++)
		{
			int iter = DOF-m;
			A_x[DOF*bc_ndnum[k]-iter] = U_new[DOF*bc_ndnum[k]-iter];
		}
    }
    // calculate r = b-A*x;
    for (int i = 1; i <= DOF*totnds; i++)
    {
        r[i] = F[i]-A_x[i];
        p[i] = 1*r[i];
    }

    int tt = 1;
    bool converged = false;
    for (tt = 1; tt <= CGiter; tt++)
    {
        U_new_iter = U_new;

        //cout << "conjugate gradient iteration " << tt << " of " << settings.galerkin_solver.CGiter << endl;

        //calculate r'*r
        double r_r = 0;
        for (int i = 1; i <= DOF*totnds; i++)
        {r_r = r_r+(r[i]*r[i]*1);}

        // calculate p'*A*p
        for(int j = 1; j <= DOF*totnds; j++)
        {A_p[j] = 0;}

        for (int s = 1; s <= tot_g_bc2; s++)
        {
            int g = g_bc2[s];
            int e = kernels[g].elem_num;

            vector<int> dof(DOF*elemnds+1);
            int kk = 1;
            for (int f = 1; f <= elemnds; f++)
            {
                for (int iter = DOF-1; iter >= 0; iter--)
                {
                    dof[kk] =  DOF*elements[e][f]-iter;
                    kk = kk + 1;
                }
            }

            for(int j = 1; j <= DOF*elemnds; j++)
            {        
		       	double product_gauss = 0; 				
                for(int i = 1; i <= DOF*elemnds; i++)
                {
                    if (DOF == 1)
                    {
                        int cm = elements[e][i];
                        int check = 0;
                        for (int h = 1; h <= bc_totnds; h++)
                        {
                            if (cm == bc_ndnum[h])
                            {check++;}
                        }
                        if (check == 0)
                        {product_gauss = product_gauss + kernels[g].MK[i][j]*p[cm];}
                    }
                    else
                    {product_gauss = product_gauss + kernels[g].MK[i][j]*p[dof[i]];}
                }
                    A_p[dof[j]] = A_p[dof[j]] + product_gauss;
            }
        }
        for (int m = 1; m <= DOF; m++)
        {
            int iter = DOF-m;
            for (int k = 1; k <= bc_totnds; k++)
            {
                A_p[DOF*bc_ndnum[k]-iter] = p[DOF*bc_ndnum[k]-iter];
            }
        }
        double p_A_p = 0;
        for (int i = 1; i <= DOF*totnds; i++)
        {p_A_p = p_A_p+(p[i]*A_p[i]);}

        double alpha=1;
        if (p_A_p == 0)
        {
            p_A_p = 1;
            cout << "WARNING: MATRIX SINGULARITY ISSUE" << endl;
        }
        else
        {alpha = r_r/p_A_p;}

        // calculate x = x+ alpha*p;
        for (int i = 1; i <= DOF*totnds; i++)
        {
            U_new[i] = U_new[i] + alpha*p[i];
            r_new[i] = r[i] - alpha*A_p[i];
        }
  
        // calculate r_new'*r_new
        double r_new_r_new = 0;
        for (int i = 1; i <= DOF*totnds; i++)
        {r_new_r_new = r_new_r_new+(r_new[i]*r_new[i]*1);}
        double betta;
        if (r_r == 0)
        {
            r_r = 1;
            cout << "ERROR: MATRIX SINGULARITY ISSUE" << endl;
            exit(0);
        }
        else
        {betta = r_new_r_new/r_r;}
        for (int i = 1; i <= DOF*totnds; i++)
        {p[i] = 1*r_new[i] + (betta*p[i]);}
        for (int i = 1; i <= DOF*totnds; i++)
        {r[i] = r_new[i];}

        double res = 0;
        for (int i = 1; i <= DOF*totnds; i++)
        {
            double diff  = abs(U_new_iter[i]-U_new[i]);
            res = res + diff;
        }
        //cout << "RES = " << res << endl;
        if (res <= conv_tol)
        {
            //cout << "Converged Solution Obtained after " << tt << " iterations..." << endl;
            converged = true;
            break;
        }
    }
    if (converged == false)
    {cout << "ERROR: No Convergence within the Specified Tolerance was Obtained for the Conjugate Gradient Solver" << endl; exit(0);}
    return U_new;
}
