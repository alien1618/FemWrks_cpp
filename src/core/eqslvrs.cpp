#include "../femwrks.h"

vector<double> solveInverse(vector<vector<double> > K, vector<double> F, int TotalPoints)
{
    vector<vector<double> > inv_K(TotalPoints+1, vector<double>(TotalPoints+1));
    vector<double> X(TotalPoints+1);
    inv_K = inv(K, TotalPoints);
    double mat_mult = 0;
    for(int i = 1; i<=TotalPoints ; i++)
    {
        for(int j = 1; j<=TotalPoints ; j++)
        {mat_mult = mat_mult + (inv_K[i][j] * F[j]);}
        X[i] = mat_mult;
        mat_mult = 0;
    }
    return X;
}
vector<double> solveTransientImplicitPCG(vector<vector<double> > K, vector<vector<double> > N, vector<double> F, vector<double> conc_now, int TotalPoints, double dt, int CGiter, double tol)
{
    cout <<"Solving system of equations..." << endl;
    vector<vector<double> > Term0(TotalPoints+1, vector<double>(TotalPoints+1,0));
    vector<vector<double> > Term1(TotalPoints+1, vector<double>(TotalPoints+1,0));
    vector<double> Term2(TotalPoints+1,0);
    vector<double> Term3(TotalPoints+1,0);
    vector<double> Term4(TotalPoints+1,0);
    vector<double> conc_future(TotalPoints+1);

    // CALCULATE: N + K*dt
    //cout <<"calculating term0" << endl;
    Term0 = multiply(K, TotalPoints, TotalPoints, dt);

    //cout <<"calculating term1" << endl;
    Term1 = add(N, Term0, TotalPoints, TotalPoints);

    // CALCULATE: F*dt
    //cout <<"calculating term2" << endl;
    Term2 = multiply(F, TotalPoints, dt);

    // CALCULATE: N * p
    //cout <<"calculating term3" << endl;
    Term3 = multiply(N, TotalPoints, TotalPoints, conc_now, TotalPoints);

    // CALCULATE: ((N * p) + (F*dt))
    //cout <<"calculating term4" << endl;
    Term4 = add(Term2, Term3, TotalPoints);

    // SOLVE THE SYSTEM OF EQUATIONS
    cout <<"Computing conjugate gradient..." << endl;
    conc_future = solvePCG(Term1, Term4, TotalPoints, CGiter, tol);
    return conc_future;
}
vector<double> solveTransientExplicitPCG(vector<vector<double> > K, vector<vector<double> > N, vector<double> F, vector<double> conc_now, int TotalPoints, double dt, int CGiter, double tol)
{
    cout <<"Solving system of equations..." << endl;
    vector<vector<double> > Term0(TotalPoints+1, vector<double>(TotalPoints+1,0));
    vector<vector<double> > Term1(TotalPoints+1, vector<double>(TotalPoints+1,0));
    vector<double> Term2(TotalPoints+1,0);
    vector<double> Term3(TotalPoints+1,0);
    vector<double> Term5(TotalPoints+1,0);
    vector<double> conc_future(TotalPoints+1);

    // CALCULATE: N + K*dt
    //cout << "calculating term0" << endl;
    Term0 = multiply(K,TotalPoints, TotalPoints, dt);

    //cout <<"calculating term1" << endl;
    Term1 = add(N, Term0, TotalPoints, TotalPoints);

    // CALCULATE: F*dt
    //cout <<"calculating term2" << endl;
    Term2 = multiply(F, TotalPoints, dt);

    // CALCULATE: (N + K*dt)*p
    //cout <<"calculating term3" << endl;
    Term3 = multiply(Term1, TotalPoints, TotalPoints, conc_now, TotalPoints);

    // CALCULATE: (N + K*dt)*p + F*dt
    //cout <<"calculating term4" << endl;
    Term5 = add(Term3, Term2, TotalPoints);

    // SOLVE THE SYSTEM OF EQUATIONS
    cout <<"Computing conjugate gradient..." << endl;
    conc_future = solvePCG( N, Term5, TotalPoints, CGiter, tol);

    return conc_future;
}

/*
--------------------------------------------------------------------------------
THIS FUNCTION SOLVES THE SYSTEM OF EQUATIONS Ax = b USING THE PRECONDITIONED CONJUGATE
GRADIENT METHOD
--------------------------------------------------------------------------------
*/
/*
vector<double> SOLVE_LinearSystemOfEquations::PCG(vector<vector<double> > A, vector<double> b, int mat_size, int CGiter, double tol)
{
    double t_start = omp_get_wtime();
    vector<double> x(mat_size+1);
    vector<double> x_old(mat_size+1);
    vector<double> A_x(mat_size+1);
    vector<double> A_p(mat_size+1);
    vector<double> r(mat_size+1);
    vector<double> r_new(mat_size+1);
    vector<double> p(mat_size+1);
    vector<vector<double> > M(mat_size+1, vector<double>(mat_size+1));
    double number;
    int threads = 1;
    #pragma omp parallel
    {threads = omp_get_num_threads();}

    //------------------------------------------
    // INITIAL GUESS FOR X
    //------------------------------------------
    x = b;

    //------------------------------------------
    // PRECONDITIONING OF MATRICES
    //------------------------------------------
#pragma omp parallel for
    for (int i = 1; i <= mat_size; i++)
    {M[i][i] = 1/A[i][i];}

#pragma omp parallel for
    for(int j = 1; j <= mat_size; j++)
    {A_x[j] = 0;}

#pragma omp parallel for
    for(int j = 1; j <= mat_size; j++)
    {
#pragma omp parallel for
        for(int i = 1; i <= mat_size; i++)
        {
            double num = A[i][j]*x[i];
            A_x[j] = A_x[j] + num;
        }
    }

#pragma omp parallel for
    for (int i = 1; i <= mat_size; i++)
    {
        r[i] = b[i]-A_x[i];
        p[i] = M[i][i]*r[i];
    }
    int t;
    double sum = 0;
    for (t = 1; t <= CGiter; t++)
    {
        x_old = x;
        // cout << "conjugate gradient iteration " << t << " of " << CGiter << endl;

        // calculate r'*r
        double r_vec[threads+1];
#pragma omp parallel
        {

            int num_threads = omp_get_num_threads();
            int this_thread = omp_get_thread_num();
            int my_start = (this_thread)*mat_size/num_threads;
            int my_end = (this_thread+1)*mat_size/num_threads;
            r_vec[this_thread] = 0;
            for (int n = my_start+1; n <= my_end; n++)
            {r_vec[this_thread] = r_vec[this_thread]+(r[n]*r[n]*M[n][n]);}
        }
        double r_r = 0;
        for(int i = 0; i <= threads; i++)
        {r_r = r_r + r_vec[i];}

        // calculate p'*A*p
#pragma omp parallel for
        for(int j = 1; j <= mat_size; j++)
        {A_p[j] = 0;}

        int i, j;
#pragma omp parallel shared(A, p, A_p) private (i,j, number)  // this line and the next are crucial. The speed up the computation X5 times
        {
#pragma omp for schedule (static)
            for(j = 1; j <= mat_size; j++)
            {
                for(i = 1; i <= mat_size; i++)
                {
                    number = A[i][j]*p[i];
                    A_p[j] = A_p[j] + number;
                }
            }
        }

        double pp[threads+1];
#pragma omp parallel
        {
            int num_threads = omp_get_num_threads();
            int this_thread = omp_get_thread_num();

            int my_start = (this_thread)*mat_size/num_threads;
            int my_end = (this_thread+1)*mat_size/num_threads;
            pp[this_thread] = 0;
            for (int n = my_start+1; n <= my_end; n++)
            {pp[this_thread] = pp[this_thread]+(p[n]*A_p[n]);}
        }
        double p_A_p = 0;
        for (int i = 0; i <= threads; i++)
        {p_A_p = p_A_p+ + pp[i];}

        double alpha=0;
        if (p_A_p == 0)
        {
            p_A_p = 1;
            cout << "WARNING: MATRIX SINGULARITY ISSUE" << endl;
        }
        else
        {alpha = r_r/p_A_p;}

        // calculate x = x+ alpha*p;
#pragma omp parallel for
        for (int i = 1; i <= mat_size; i++)
        {
            x[i] = x[i] + alpha*p[i];
            r_new[i] = r[i] - alpha*A_p[i];
        }

        double rr[threads+1];
#pragma omp parallel
        {
            int num_threads = omp_get_num_threads();
            int this_thread = omp_get_thread_num();

            int my_start = (this_thread)*mat_size/num_threads;
            int my_end = (this_thread+1)*mat_size/num_threads;
            rr[this_thread] = 0;
            for (int n = my_start+1; n <= my_end; n++)
            {rr[this_thread] = rr[this_thread]+(r_new[n]*r_new[n]*M[n][n]);}
        }
        double r_new_r_new = 0;
        for (int i = 0; i <= threads; i++)
        {r_new_r_new = r_new_r_new + rr[i];}


        double betta;
        if (r_r == 0)
        {
            r_r = 1;
            cout << "ERROR: MATRIX SINGULARITY ISSUE" << endl;
            exit(0);
        }
        else
        {betta = r_new_r_new/r_r;}

#pragma omp parallel for
        for (int i = 1; i <= mat_size; i++)
        {p[i] = M[i][i]*r_new[i] + (betta*p[i]);}

#pragma omp parallel for
        for (int i = 1; i <= mat_size; i++)
        {r[i] = r_new[i];}

        double s[threads+1];
#pragma omp parallel
        {
            int num_threads = omp_get_num_threads();
            int this_thread = omp_get_thread_num();

            int my_start = (this_thread)*mat_size/num_threads;
            int my_end = (this_thread+1)*mat_size/num_threads;
            s[this_thread] = 0;
            for (int n = my_start+1; n <= my_end; n++)
            {s[this_thread] = s[this_thread]+abs(x_old[n]-x[n]);}
        }
        sum = 0;
        for (int i = 0; i <= threads; i++)
        {sum = sum + s[i];}

        if (sum <= tol)
        {break;}
    }
    if (t < CGiter)
    {
        cout << "Converged Matrix Solution Obtained after " << t << " iterations..." << endl;
        converged = true;
    }
    else
    {
        cout << "No Convergence within the Specified Tolerance was Obtained for the Conjugate Gradient Solver" << endl;
        converged = false;
        exit(0);
    }
    double t_end = omp_get_wtime() - t_start;
    cout << "Total time elapsed for computing conjugate gradient = " << t_end << endl;
    return x;

}
*/
vector<double> solvePCG(vector<vector<double> > A, vector<double> b, int mat_size, int CGiter, double tol)
{
    //Serial Code with no OpenMP
    vector<double> x(mat_size+1);
    vector<double> x_old(mat_size+1);
    vector<double> A_x(mat_size+1);
    vector<double> A_p(mat_size+1);
    vector<double> r(mat_size+1);
    vector<double> r_new(mat_size+1);
    vector<double> p(mat_size+1);
    double number;

    // initial guess of x
    x = b;

    // calculate r = b-A*x;
    for(int j = 1; j <= mat_size; j++)
    {A_x[j] = 0;}
    for(int j = 1; j <= mat_size; j++)
    {
        for(int i = 1; i <= mat_size; i++)
        {
            number = A[i][j]*x[i];
            A_x[j] = A_x[j] + number;
        }
    }
    for (int i = 1; i <= mat_size; i++)
    {
        r[i] = b[i]-A_x[i];
        p[i] = (1/A[i][i])*r[i];
    }
    int t;
    for (t = 1; t <= CGiter; t++)
    {

        x_old = x;

        //cout << "conjugate gradient iteration " << t << " of " << CGiter << endl;

        //calculate r'*r
        double r_r = 0;
        for (int i = 1; i <= mat_size; i++)
        {r_r = r_r+(r[i]*r[i]*(1/A[i][i]));}

        // calculate p'*A*p
        for(int j = 1; j <= mat_size; j++)
        {A_p[j] = 0;}
        for(int j = 1; j <= mat_size; j++)
        {
            for(int i = 1; i <= mat_size; i++)
            {
                number = A[i][j]*p[i];
                A_p[j] = A_p[j] + number;
            }
        }

        double p_A_p = 0;
        for (int i = 1; i <= mat_size; i++)
        {p_A_p = p_A_p+(p[i]*A_p[i]);}

        double alpha=0;
        if (p_A_p == 0)
        {
            p_A_p = 1;
            cout << "WARNING: MATRIX SINGULARITY ISSUE" << endl;
        }
        else
        {alpha = r_r/p_A_p;}

        // calculate x = x+ alpha*p;
        for (int i = 1; i <= mat_size; i++)
        {
            x[i] = x[i] + alpha*p[i];
            r_new[i] = r[i] - alpha*A_p[i];
        }
        //if (r_new[1] <= 0.00001)
        //{ break;}

        // calculate r_new'*r_new
        double r_new_r_new = 0;
        for (int i = 1; i <= mat_size; i++)
        {r_new_r_new = r_new_r_new+(r_new[i]*r_new[i]*(1/A[i][i]));}
        double betta;
        if (r_r == 0)
        {
            r_r = 1;
            cout << "ERROR: MATRIX SINGULARITY ISSUE" << endl;
            exit(0);
        }
        else
        {betta = r_new_r_new/r_r;}
        for (int i = 1; i <= mat_size; i++)
        {p[i] = (1/A[i][i])*r_new[i] + (betta*p[i]);}
        for (int i = 1; i <= mat_size; i++)
        {r[i] = r_new[i];}

        double sum = 0;
        for (int i = 1; i <= mat_size; i++)
        {
            double diff  = abs(x_old[i]-x[i]);
            sum = sum + diff;
        }
        if (sum <= tol)
        {break;}
    }
    if (t < CGiter)
    {cout << "Converged Matrix Solution Obtained after " << t << " iterations..." << endl;}
    else
    {cout << "ERROR: No Convergence within the Specified Tolerance was Obtained for the Conjugate Gradient Solver" << endl; exit(0);}
    return x;

}
