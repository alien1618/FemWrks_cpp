#include "../femwrks.h"

double mod(double number, double divnum)
{
    double ans = number/divnum;
    double remainder = number - double(int(ans) * divnum);
    return remainder;
}

double maxnum(double a, double b)
{
    double answer=0;
    if (a > b)
    {answer = a;}
    else if (b > a)
    {answer = b;}
    else if (b == a)
    {answer = b;}
    return answer;
}
double minnum(double a, double b)
{
    double answer=0;
    if (a < b)
    {answer = a;}
    else if (b < a)
    {answer = b;}
    else if (b == a)
    {answer = b;}
    return answer;
}

double random0to1()
{return rand()/ double(RAND_MAX);}
void randomSeed()
{srand(time(0));}


double det_2X2(vector<vector<double> > mat)
{
    double det;
    det = (mat[1][1] * mat[2][2]) - (mat[2][1]* mat[1][2]);
    return det;
}

double det_3X3(vector<vector<double> > mat)
{
    double det;
    det = mat[1][1] * ((mat[2][2]*mat[3][3])-(mat[3][2]*mat[2][3]));
    det = det - (mat[2][1] * ((mat[1][2]*mat[3][3])-(mat[3][2]*mat[1][3])));
    det = det + (mat[3][1] * ((mat[1][2]*mat[2][3])-(mat[2][2]*mat[1][3])));
    return det;
}
vector<vector<double> > inverse_2X2(vector<vector<double> > A)
{
    vector<vector<double> > A_inv(3, vector<double> (3));
    double num = 1/det_2X2(A);
    A_inv[1][1] = num*A[2][2];
    A_inv[2][1] = -num*A[2][1];
    A_inv[1][2] = -num*A[1][2];
    A_inv[2][2] = num*A[1][1];
    return A_inv;
}
vector<vector<double> > inverse_3X3(vector<vector<double> > A)
{
    vector<vector<double> > A_inv(4, vector<double> (4));
    vector<vector<double> > T(4, vector<double> (4));
    double num = 1/det_3X3(A);

    for (int j = 1; j<= 3; j++)
    {
        for (int i = 1; i <= 3; i++)
        {T[i][j] = A[j][i];}
    }

    double T11 = T[2][2]*T[3][3] - T[3][2]*T[2][3];
    double T21 = -(T[1][2]*T[3][3] - T[3][2]*T[1][3]);
    double T31 = T[1][2]*T[2][3] - T[2][2]*T[1][3];
    double T12 = -(T[2][1]*T[3][3] - T[3][1]*T[2][3]);
    double T22 = T[1][1]*T[3][3] - T[3][1]*T[1][3];
    double T32 = -(T[1][1]*T[2][3] - T[2][1]*T[1][3]);
    double T13 = T[2][1]*T[3][2] - T[3][1]*T[2][2];
    double T23 = -(T[1][1]*T[3][2] - T[3][1]*T[1][2]);
    double T33 = T[1][1]*T[2][2] - T[2][1]*T[1][2];

    A_inv[1][1] = num*T11;
    A_inv[2][1] = num*T21;
    A_inv[3][1] = num*T31;
    A_inv[1][2] = num*T12;
    A_inv[2][2] = num*T22;
    A_inv[3][2] = num*T32;
    A_inv[1][3] = num*T13;
    A_inv[2][3] = num*T23;
    A_inv[3][3] = num*T33;
    return A_inv;
}

// Matrix A + Matrix B
vector<double> add(vector<double> A, vector<double> B, int size)
{
    vector<double> C(size+1);
    for (int i = 1; i <= size; i++)
    {C[i] = A[i] + B[i];}
    return C;
}

// Matrix A + Matrix B
vector<vector<double> > add(vector<vector<double> > A, vector<vector<double> > B, int Xclmns, int rows)
{
    vector<vector<double> > C(Xclmns+1, vector<double>(rows+1));
#pragma omp parallel for
    for (int j = 1; j <= rows; j++)
    {
#pragma omp parallel for
        for (int i = 1; i <= Xclmns; i++)
        {C[i][j] = A[i][j] + B[i][j];}
    }
    return C;
}
// Matrix A * Matrix B
vector<vector<double> > multiply(vector<vector<double> > A, int A_clmns, int A_rows, vector<vector<double> > B, int B_clmns, int B_rows)
{
    vector<vector<double> > C(B_clmns+1, vector<double>(A_rows+1));
    if (A_clmns != B_rows)
    {cout << endl << "ERROR. matrix multiplication failure. matrix size mismatch" << endl;exit(0);}
    else
    {
        double sum=0;
        for (int j = 1; j <= A_rows; j++)
        {
            for (int k = 1; k <= B_clmns; k++)
            {
                sum = 0;
                for (int i = 1; i <= A_clmns; i++)
                {sum = sum + A[i][j] * B[k][i];}
                C[k][j] = sum;
            }
        }
    }
    return C;
}

// number * Matrix A
vector<vector<double> > multiply(vector<vector<double> > A, int A_Xclmns, int A_rows, double number)
{
#pragma omp parallel for
    for (int j = 1; j <= A_rows; j++)
    {
#pragma omp parallel for
        for (int i = 1; i <= A_Xclmns; i++)
        {
            A[i][j] = A[i][j] * number;
        }
    }
    return A;
}

// number * Matrix A
vector<double> multiply(vector<double> A, int A_size, double number)
{
#pragma omp parallel for
    for (int i = 1; i <= A_size; i++)
    {A[i] = A[i] * number;}
    return A;
}

// Matrix A' * Matrix A
vector<vector<double> > multiply_transp(vector<double> A, vector<double> B, int size)
{
    vector<vector<double> >  C(size+1, vector<double>(size+1));
    int i,j;
#pragma omp parallel shared(A,B,C) private(i,j)
    {
#pragma omp for schedule (static)
        for (j = 1; j <= size; j++)
        {
            for (i = 1; i <= size; i++)
            {C[i][j] = A[j] * B[i];}
        }
    }
    return C;
}
vector<double> multiply(vector<vector<double> > A, int A_Xclmns, int A_Yrows, vector<double> B, int B_size)
{
    vector<double> C(A_Yrows+1);
    if (A_Xclmns != B_size)
    {cout << endl << "ERROR. matrix multiplication failure. matrix size mismatch" << endl;exit(0);}
    else
    {

#pragma omp parallel shared(A,B,C)
        {
            double mat_mult = 0;
#pragma omp for schedule (static)
            for(int j = 1; j<= A_Yrows ; j++)
            {
                for(int i = 1; i<= B_size ; i++)
                {mat_mult = mat_mult + (A[i][j] * B[i]);}
                C[j] = mat_mult;
                mat_mult = 0;
            }
        }
    }
    return C;
}
// transpose of Matrix A
vector<vector<double> > transpose(vector<vector<double> > A, int A_Xclmns, int A_rows)
{
    vector<vector<double> > C(A_rows+1, vector<double>(A_Xclmns+1));
#pragma omp parallel for
    for (int j = 1; j <= A_rows; j++)
    {
        for (int i = 1; i <= A_Xclmns; i++)
        {C[j][i] = A[i][j];}
    }
    return C;
}
vector<vector<double> > inv(vector<vector<double> > a, int size)
{
    vector<vector<double> > b(size+1, vector<double>(size+1));
    double s,t;
#pragma omp parallel for
    for (int j = 1; j <= size; j++)
    {
#pragma omp parallel for
        for(int i = 1; i <= size; i++)
        {b[i][j] = 0;}
    }
#pragma omp parallel for
    for (int i = 1; i<= size; i++)
    {b[i][i] = 1;}
    //The following code actually performs the matrix inversion
    for (int j = 1; j <= size; j++)
    {
        for (int i = j; i<= size; i++)
        {
            if (a[i][j] != 0)
            {
                for (int k = 1; k<= size; k++)
                {
                    s = a[j][k];
                    a[j][k] = a[i][k];
                    a[i][k] = s;
                    s = b[j][k];
                    b[j][k] = b[i][k];
                    b[i][k] = s;
                }
                t = 1/a[j][j];

                for (int k = 1; k <= size; k++)
                {
                    a[j][k] = t * a[j][k];
                    b[j][k] = t * b[j][k];
                }
                for (int L = 1; L <= size; L++)
                {
                    if (L != j)
                    {
                        t = -a[L][j];
                        for (int k = 1; k<= size; k++)
                        {
                            a[L][k] = a[L][k] + t * a[j][k];
                            b[L][k] = b[L][k] + t * b[j][k];
                        }
                    }
                }
            }
            break;
        }
    }
    a.clear();
    return b;
}


