#include "../femwrks.h"

vector<KERNEL> computeInterpolants(int dim, vector<POINT> points, vector<vector<int> > elements, int totelems, int elemnds, GAUSS gauss)
{
    int c= 1;
    int m = 1;
    int TotalPoints = 0;
    vector<KERNEL> kernels;
    TotalPoints = gauss.element_integpoints*totelems;
    kernels.resize(TotalPoints+1);
    for (int k = 1; k <= totelems; k++)
    {
        for (int i = 1; i <= gauss.element_integpoints; i++)
        {
            kernels[m].computeLagrange(dim, points, elements[k], elemnds, gauss.eps[i],gauss.eta[i], gauss.zta[i]);
            kernels[m].quad_weight = gauss.w[i]*abs(kernels[m].det_jac);
            kernels[m].elem_num = k;
            
            kernels[m].dNdNx = multiply_transp(kernels[m].dNdx, kernels[m].dNdx, elemnds);
            kernels[m].dNdNy = multiply_transp(kernels[m].dNdy, kernels[m].dNdy, elemnds);
            kernels[m].dNdNz = multiply_transp(kernels[m].dNdz, kernels[m].dNdz, elemnds);
            kernels[m].NdNx = multiply_transp(kernels[m].N, kernels[m].dNdx, elemnds);
            kernels[m].NdNy = multiply_transp(kernels[m].N, kernels[m].dNdy, elemnds);
            kernels[m].NdNz = multiply_transp(kernels[m].N, kernels[m].dNdz, elemnds);
            kernels[m].NN = multiply_transp(kernels[m].N, kernels[m].N, elemnds);

            m++;
        }
        if (k/1000 == c)
        {
            c++;
            cout << "Generated quad points inside " << k << " of " << totelems << " surfaces..." << endl;
        }
    }        
    return kernels;
}

void KERNEL::computeLagrange(int dim, vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps, double eta, double zta)
{
    switch (dim)
    {
        case 1:
        {lineO1(points, nbrnds, TotNbrNds, eps);}break;
        case 2:
        {
            switch (TotNbrNds)
            {
                case 4:
                    Q4(points, nbrnds, TotNbrNds,eps,eta);break;
                case 8:
                    Q8(points,nbrnds, TotNbrNds,eps,eta);break;
                case 3:
                    T3(points,nbrnds, TotNbrNds,eps,eta);break;
                case 6:
                    T6(points,nbrnds, TotNbrNds,eps,eta);break;
            }
        }break;
        case 3:
        {
            switch (TotNbrNds)
            {
                case 8:
                    H8(points,nbrnds, TotNbrNds,eps, eta, zta);break;
                case 16:
                    H20(points,nbrnds, TotNbrNds,eps,eta, zta);break;
                case 4:
                    T4(points,nbrnds, TotNbrNds,eps,eta, zta);break;
                case 10:
                    T10(points,nbrnds, TotNbrNds,eps,eta, zta);break;
            }
        }break;
    }
}
void KERNEL::lineO1(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps)
{
    if (TotNbrNds != 2)
    {cout << "ERROR in KERNEL::lineO1. Total support points incompatible with element type" << endl; exit(0);}

    double x1 = points[nbrnds[1]].x;
    double y1 = points[nbrnds[1]].y;

    double x2 = points[nbrnds[2]].x;
    double y2 = points[nbrnds[2]].y;

    double ycen = (y2+y1)/2;
    double xcen = (x2+x1)/2;
    double jcob_y = abs((y2-y1)/2);
    double jcob_x = abs((x2-x1)/2);

    //------------------------------------------
    // generating interpolation functions
    //------------------------------------------
    double x = xcen-eps*jcob_x;
    double y = ycen-eps*jcob_y;

    vector<double> N(3);
    vector<double> dNdx(3);
    vector<double> dNdy(3);
    vector<double> dNdz(3);
    if (y1 - y2 == 0)
    {
        N[1] = (x-x2)/(x1-x2);
        dNdx[1] = 1;
        dNdy[1] = 0;
        dNdz[1] = 0;
        N[2] = 1-N[1];
        dNdx[2] = -1;
        dNdy[2] = 0;
        dNdz[2] = 0;
    }
    else
    {
        N[1] = (y-y2)/(y1-y2);
        dNdy[1] = 1;
        dNdx[1] = 0;
        N[2] = 1-N[1];
        dNdy[2] = -1;
        dNdx[2] = 0;
        dNdz[1] = 0;
        dNdz[2] = 0;
    }
}
void KERNEL::H8(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps,double eta, double zta)
{
    if (TotNbrNds != 8)
    {cout << "ERROR in KERNEL::H8. Total support points incompatible with element type" << endl; exit(0);}

    N.resize(TotNbrNds+1);
    vector<double> dNdeps(TotNbrNds+1);
    vector<double> dNdeta(TotNbrNds+1);
    vector<double> dNdzta(TotNbrNds+1);

    N[1] = 0.125*(1-eps)*(1-eta)*(1+zta);
    N[2] = 0.125*(1+eps)*(1-eta)*(1+zta);
    N[3] = 0.125*(1+eps)*(1+eta)*(1+zta);
    N[4] = 0.125*(1-eps)*(1+eta)*(1+zta);
    N[5] = 0.125*(1-eps)*(1-eta)*(1-zta);
    N[6] = 0.125*(1+eps)*(1-eta)*(1-zta);
    N[7] = 0.125*(1+eps)*(1+eta)*(1-zta);
    N[8] = 0.125*(1-eps)*(1+eta)*(1-zta);

    dNdeps[1] = -0.125*(1-eta)*(1+zta);
    dNdeta[1] = -0.125*(1-eps)*(1+zta);
    dNdzta[1] = 0.125*(1-eps)*(1-eta);

    dNdeps[2] = 0.125*(1-eta)*(1+zta);
    dNdeta[2] = -0.125*(1+eps)*(1+zta);
    dNdzta[2] = 0.125*(1+eps)*(1-eta);

    dNdeps[3] = 0.125*(1+eta)*(1+zta);
    dNdeta[3] = 0.125*(1+eps)*(1+zta);
    dNdzta[3] = 0.125*(1+eps)*(1+eta);

    dNdeps[4] = -0.125*(1+eta)*(1+zta);
    dNdeta[4] = 0.125*(1-eps)*(1+zta);
    dNdzta[4] = 0.125*(1-eps)*(1+eta);

    dNdeps[5] = -0.125*(1-eta)*(1-zta);
    dNdeta[5] = -0.125*(1-eps)*(1-zta);
    dNdzta[5] = -0.125*(1-eps)*(1-eta);

    dNdeps[6] = 0.125*(1-eta)*(1-zta);
    dNdeta[6] = -0.125*(1+eps)*(1-zta);
    dNdzta[6] = -0.125*(1+eps)*(1-eta);

    dNdeps[7] = 0.125*(1+eta)*(1-zta);
    dNdeta[7] = 0.125*(1+eps)*(1-zta);
    dNdzta[7] = -0.125*(1+eps)*(1+eta);

    dNdeps[8] = -0.125*(1+eta)*(1-zta);
    dNdeta[8] = 0.125*(1-eps)*(1-zta);
    dNdzta[8] = -0.125*(1-eps)*(1+eta);

    jacobianVolume(points, nbrnds, TotNbrNds, dNdeps, dNdeta, dNdzta);
}
void KERNEL::H20(vector<POINT> points, vector<int> nbrnds, int TotNbrNds,double eps,double eta, double zta)
{
    if (TotNbrNds != 20)
    {cout << "ERROR in KERNEL::H20. Total support points incompatible with element type" << endl; exit(0);}
    N.resize(TotNbrNds+1);
    vector<double> dNdeps(TotNbrNds+1);
    vector<double> dNdeta(TotNbrNds+1);
    vector<double> dNdzta(TotNbrNds+1);

    N[1] = 0.125*(1-eps)*(1-eta)*(1-zta)*(-eps-eta-zta-2);
    N[2] = 0.125*(1-eps)*(1-eta)*(1-zta)*(eps-eta-zta-2);
    N[3] = 0.125*(1+eps)*(1+eta)*(1-zta)*(eps+eta-zta-2);
    N[4] = 0.125*(1-eps)*(1+eta)*(1-zta)*(-eps+eta-zta-2);
    N[5] = 0.125*(1-eps)*(1-eta)*(1+zta)*(-eps-eta+zta-2);
    N[6] = 0.125*(1+eps)*(1-eta)*(1+zta)*(eps-eta+zta-2);
    N[7] = 0.125*(1+eps)*(1+eta)*(1+zta)*(eps+eta+zta-2);
    N[8] = 0.125*(1-eps)*(1+eta)*(1+zta)*(-eps+eta+zta-2);
    N[9] = 0.125*2*(1-eps*eps)*(1-eta)*(1-zta);
    N[10] = 0.125*2*(1+eps)*(1-eta*eta)*(1-zta);
    N[11] = 0.125*2*(1-eps*eps)*(1+eta)*(1-zta);
    N[12] = 0.125*2*(1-eps)*(1-eta*eta)*(1-zta);
    N[13] = 0.125*2*(1-eps)*(1-eta)*(1-zta*zta);
    N[14] = 0.125*2*(1+eps)*(1-eta)*(1-zta*zta);
    N[15] = 0.125*2*(1+eps)*(1+eta)*(1-zta*zta);
    N[16] = 0.125*2*(1-eps)*(1+eta)*(1-zta*zta);
    N[17] = 0.125*2*(1-eps*eps)*(1-eta)*(1+zta);
    N[18] = 0.125*2*(1+eps)*(1-eta*eta)*(1+zta);
    N[19] = 0.125*2*(1-eps*eps)*(1+eta)*(1+zta);
    N[20] = 0.125*2*(1-eps)*(1-eta*eta)*(1+zta);

    dNdeps[1] = -0.125*(1-eta)*(1-zta)*(-eps-eta-zta-2)-0.125*(1-eps)*(1-eta)*(1-zta);
    dNdeta[1] = -0.125*(1-eps)*(1-zta)*(-eps-eta-zta-2)-0.125*(1-eps)*(1-eta)*(1-zta);
    dNdzta[1] = -0.125*(1-eps)*(1-eta)*(-eps-eta-zta-2)-0.125*(1-eps)*(1-eta)*(1-zta);

    dNdeps[2] = -0.125*(1-eta)*(1-zta)*(eps-eta-zta-2)+0.125*(1-eps)*(1-eta)*(1-zta);
    dNdeta[2] = -0.125*(1-eps)*(1-zta)*(eps-eta-zta-2)-0.125*(1-eps)*(1-eta)*(1-zta);
    dNdzta[2] = -0.125*(1-eps)*(1-eta)*(eps-eta-zta-2)-0.125*(1-eps)*(1-eta)*(1-zta);

    dNdeps[3] = 0.125*(1+eta)*(1-zta)*(eps+eta-zta-2)+0.125*(1+eps)*(1+eta)*(1-zta);
    dNdeta[3] = 0.125*(1+eps)*(1-zta)*(eps+eta-zta-2)+0.125*(1+eps)*(1+eta)*(1-zta);
    dNdzta[3] = -0.125*(1+eps)*(1+eta)*(eps+eta-zta-2)-0.125*(1+eps)*(1+eta)*(1-zta);

    dNdeps[4] = -0.125*(1+eta)*(1-zta)*(-eps+eta-zta-2)-0.125*(1-eps)*(1+eta)*(1-zta);
    dNdeta[4] =  0.125*(1-eps)*(1-zta)*(-eps+eta-zta-2)+ 0.125*(1-eps)*(1+eta)*(1-zta);
    dNdzta[4] =  -0.125*(1-eps)*(1+eta)*(-eps+eta-zta-2)-0.125*(1-eps)*(1+eta)*(1-zta);

    dNdeps[5] = -0.125*(1-eta)*(1+zta)*(-eps-eta+zta-2)-0.125*(1-eps)*(1-eta)*(1+zta);
    dNdeta[5] = -0.125*(1-eps)*(1+zta)*(-eps-eta+zta-2)-0.125*(1-eps)*(1-eta)*(1+zta);
    dNdzta[5] = 0.125*(1-eps)*(1-eta)*(-eps-eta+zta-2)+0.125*(1-eps)*(1-eta)*(1+zta);

    dNdeps[6] = 0.125*(1-eta)*(1+zta)*(eps-eta+zta-2)+0.125*(1+eps)*(1-eta)*(1+zta);
    dNdeta[6] = -0.125*(1+eps)*(1+zta)*(eps-eta+zta-2)-0.125*(1+eps)*(1-eta)*(1+zta);
    dNdzta[6] = 0.125*(1+eps)*(1-eta)*(eps-eta+zta-2)+0.125*(1+eps)*(1-eta)*(1+zta);

    dNdeps[7] = 0.125*(1+eta)*(1+zta)*(eps+eta+zta-2)+0.125*(1+eps)*(1+eta)*(1+zta);
    dNdeta[7] = 0.125*(1+eps)*(1+zta)*(eps+eta+zta-2)+0.125*(1+eps)*(1+eta)*(1+zta);
    dNdzta[7] = 0.125*(1+eps)*(1+eta)*(eps+eta+zta-2)+0.125*(1+eps)*(1+eta)*(1+zta);

    dNdeps[8] = -0.125*(1+eta)*(1+zta)*(-eps+eta+zta-2)-0.125*(1-eps)*(1+eta)*(1+zta);
    dNdeta[8] = 0.125*(1-eps)*(1+zta)*(-eps+eta+zta-2)+0.125*(1-eps)*(1+eta)*(1+zta);
    dNdzta[8] = 0.125*(1-eps)*(1+eta)*(-eps+eta+zta-2)+0.125*(1-eps)*(1+eta)*(1+zta);

    dNdeps[9] = 0.125*2*(-2*eps)*(1-eta)*(1-zta);
    dNdeta[9] = -0.125*2*(1-eps*eps)*(1-zta);
    dNdzta[9] = -0.125*2*(1-eps*eps)*(1-eta);

    dNdeps[10] = 0.125*2*(1-eta*eta)*(1-zta);
    dNdeta[10] = 0.125*2*(1+eps)*(-2*eta)*(1-zta);
    dNdzta[10] = -0.125*2*(1+eps)*(1-eta*eta);

    dNdeps[11] = 0.125*2*(-2*eps)*(1+eta)*(1-zta);
    dNdeta[11] = 0.125*2*(1-eps*eps)*(1-zta);
    dNdzta[11] = -0.125*2*(1-eps*eps)*(1+eta);

    dNdeps[12] =  -0.125*2*(1-eta*eta)*(1-zta);
    dNdeta[12] =  0.125*2*(1-eps)*(-2*eta)*(1-zta);
    dNdzta[12] =  -0.125*2*(1-eps)*(1-eta*eta);

    dNdeps[13] = -0.125*2*(1-eta)*(1-zta*zta);
    dNdeta[13] = -0.125*2*(1-eps)*(1-zta*zta);
    dNdzta[13] = 0.125*2*(1-eps)*(1-eta)*(-2*zta);

    dNdeps[14] = 0.125*2*(1-eta)*(1-zta*zta);
    dNdeta[14] =  -0.125*2*(1+eps)*(1-zta*zta);
    dNdzta[14] =  0.125*2*(1+eps)*(1-eta)*(-2*zta);

    dNdeps[15] = 0.125*2*(1+eta)*(1-zta*zta);
    dNdeta[15] = 0.125*2*(1+eps)*(1-zta*zta);
    dNdzta[15] = 0.125*2*(1+eps)*(1+eta)*(-2*zta);

    dNdeps[16] = -0.125*2*(1+eta)*(1-zta*zta);
    dNdeta[16] = 0.125*2*(1-eps)*(1-zta*zta);
    dNdzta[16] = 0.125*2*(1-eps)*(1+eta)*(-2*zta);

    dNdeps[17] = 0.125*2*(-2*eps)*(1-eta)*(1+zta);
    dNdeta[17] = -0.125*2*(1-eps*eps)*(1+zta);
    dNdzta[17] = 0.125*2*(1-eps*eps)*(1-eta);

    dNdeps[18] = 0.125*2*(1-eta*eta)*(1+zta);
    dNdeta[18] = 0.125*2*(1+eps)*(-2*eta)*(1+zta);
    dNdzta[18] = 0.125*2*(1+eps)*(1-eta*eta);

    dNdeps[19] = 0.125*2*(-2*eps)*(1+eta)*(1+zta);
    dNdeta[19] = 0.125*2*(1-eps*eps)*(1+zta);
    dNdzta[19] = 0.125*2*(1-eps*eps)*(1+eta);

    dNdeps[20] = -0.125*2*(1-eta*eta)*(1+zta);
    dNdeta[20] = 0.125*2*(1-eps)*(-2*eta)*(1+zta);
    dNdzta[20] = 0.125*2*(1-eps)*(1-eta*eta);

    jacobianVolume(points, nbrnds, TotNbrNds, dNdeps, dNdeta, dNdzta);
}
void  KERNEL::T4(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps,double eta, double zta)
{
    if (TotNbrNds != 4)
    {cout << "ERROR in KERNEL::T4. Total support points incompatible with element type" << endl; exit(0);}

    N.resize(TotNbrNds+1);
    vector<double> dNdeps(TotNbrNds+1);
    vector<double> dNdeta(TotNbrNds+1);
    vector<double> dNdzta(TotNbrNds+1);

    N[1] = 1-eps-eta-zta;
    N[2] = eps;
    N[3] = eta;
    N[4] = zta;

    dNdeps[1] = -1.0;
    dNdeta[1] = -1.0;
    dNdzta[1] = -1.0;

    dNdeps[2] = 1.0;
    dNdeta[2] = 0.0;
    dNdzta[2] = 0.0;

    dNdeps[3] = 0.0;
    dNdeta[3] = 1.0;
    dNdzta[3] = 0.0;

    dNdeps[4] = 0.0;
    dNdeta[4] = 0.0;
    dNdzta[4] = 1.0;

    jacobianVolume(points, nbrnds, TotNbrNds, dNdeps, dNdeta, dNdzta);

}
void  KERNEL::T10(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps,double eta, double zta)
{
    if (TotNbrNds != 10)
    {cout << "ERROR in KERNEL::T10. Total support points incompatible with element type" << endl; exit(0);}

    N.resize(TotNbrNds+1);
    vector<double> dNdeps(TotNbrNds+1);
    vector<double> dNdeta(TotNbrNds+1);
    vector<double> dNdzta(TotNbrNds+1);

    N[1] = (1 - eps - eta - zta)*(2*(1 - eps - eta - zta)-1);
    N[2] = eps*(2*eps-1);
    N[3] = eta*(2*eta-1);
    N[4] = zta*(2*zta-1);
    N[5] = 4*(1 - eps - eta - zta)*eps;
    N[6] = 4*eps*eta;
    N[7] = 4*(1 - eps - eta - zta)*eta;
    N[8] = 4*(1 - eps - eta - zta)*zta;
    N[9] = 4*eps*zta;
    N[10] = 4*eta*zta;

    dNdeps[1] = -1*(2*(1 - eps - eta - zta)-1)+(1 - eps - eta - zta)*(-2);
    dNdeta[1] = -1*(2*(1 - eps - eta - zta)-1)+(1 - eps - eta - zta)*(-2);
    dNdzta[1] = -1*(2*(1 - eps - eta - zta)-1)+(1 - eps - eta - zta)*(-2);

    dNdeps[2] = 4*eps-1;
    dNdeta[2] = 0;
    dNdzta[2] = 0;

    dNdeps[3] = 0;
    dNdeta[3] = 4*eta-1;
    dNdzta[3] = 0;

    dNdeps[4] = 0;
    dNdeta[4] = 0;
    dNdzta[4] = 4*zta-1;

    dNdeps[5] = -4*eps+4*(1 - eps - eta - zta);
    dNdeta[5] = -4*eps;
    dNdzta[5] = -4*eps;

    dNdeps[6] = 4*eta;
    dNdeta[6] = 4*eps;
    dNdzta[6] = 0;

    dNdeps[7] = -4*eta;
    dNdeta[7] = -4*eta+4*(1 - eps - eta - zta);
    dNdzta[7] = -4*eta;

    dNdeps[8] = -4*zta;
    dNdeta[8] = -4*zta;
    dNdzta[8] = -4*zta+4*(1 - eps - eta - zta);

    dNdeps[9] = 4*zta;
    dNdeta[9] = 0;
    dNdzta[9] = 4*eps;

    dNdeps[10] = 0;
    dNdeta[10] = 4*zta;
    dNdzta[10] = 4*eta;

    jacobianVolume(points, nbrnds, TotNbrNds, dNdeps, dNdeta, dNdzta);
}
void KERNEL::Q4(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps,double eta)
{
    if (TotNbrNds != 4)
    {cout << "ERROR in KERNEL::Q4. Total support points incompatible with element type" << endl; exit(0);}

    N.resize(TotNbrNds+1);
    vector<double> dNdeps(TotNbrNds+1);
    vector<double> dNdeta(TotNbrNds+1);

    N[1] = 0.25*(1-eps)*(1-eta);
    N[2] = 0.25*(1+eps)*(1-eta);
    N[3] = 0.25*(1+eps)*(1+eta);
    N[4] = 0.25*(1-eps)*(1+eta);

    dNdeps[1] = -0.25*(1-eta);
    dNdeta[1] = -0.25*(1-eps);
    
    dNdeps[2] = 0.25*(1-eta);
    dNdeta[2] = -0.25*(1+eps);
    
    dNdeps[3] = 0.25*(1+eta);
    dNdeta[3] = 0.25*(1+eps);
    
    dNdeps[4] = -0.25*(1+eta);
    dNdeta[4] = 0.25*(1-eps);

    jacobianSurface(points, nbrnds, TotNbrNds, dNdeps, dNdeta);
}
void KERNEL::Q8(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps,double eta)
{
    if (TotNbrNds != 8)
    {cout << "ERROR in KERNEL::Q8. Total support points incompatible with element type" << endl; exit(0);}

    N.resize(TotNbrNds+1);
    vector<double> dNdeps(TotNbrNds+1);
    vector<double> dNdeta(TotNbrNds+1);

    N[1] = 0.25*(eps-1)*(1-eta)*(eps+eta+1);
    N[2] = 0.25*(eps+1)*(1-eta)*(eps-eta-1);
    N[3] = 0.25*(1+eps)*(1+eta)*(eps+eta-1);
    N[4] = 0.25*(1-eps)*(1+eta)*(-eps+eta-1);
    N[5] = 0.5*(1-eps*eps)*(1-eta);
    N[6] = 0.5*(1+eps)*(1-eta*eta);
    N[7] = 0.5*(1-eps*eps)*(1+eta);
    N[8] = 0.5*(1-eps)*(1-eta*eta);

    dNdeps[1] = 0.25*(eps-1)*(1-eta) + 0.25*(1-eta)*(eps+eta+1);
    dNdeta[1] = -0.25*(eps-1)*(eps+eta+1) + 0.25*(eps-1)*(1-eta);

    dNdeps[2] = 0.25*(1-eta)*(eps-eta-1) + 0.25*(eps+1)*(1-eta);
    dNdeta[2] = -0.25*(eps+1)*(eps-eta-1) - 0.25*(eps+1)*(1-eta);

    dNdeps[3] = 0.25*(1+eta)*(eps+eta-1) + 0.25*(1+eps)*(1+eta);
    dNdeta[3] = 0.25*(1+eps)*(eps+eta-1) + 0.25*(1+eps)*(1+eta);

    dNdeps[4] = -0.25*(1+eta)*(-eps+eta-1) - 0.25*(1-eps)*(1+eta);
    dNdeta[4] = 0.25*(1-eps)*(-eps+eta-1) + 0.25*(1-eps)*(1+eta);

    dNdeps[5] = -0.5*(1-eta)*(2*eps);
    dNdeta[5] = -0.5*(1-eps*eps);

    dNdeps[6] = 0.5*(1-eta*eta);
    dNdeta[6] = -0.5*(1+eps)*(2*eta);

    dNdeps[7]= -0.5*(1+eta)*(2*eps);
    dNdeta[7]= 0.5*(1-eps*eps);

    dNdeps[8] = -0.5*(1-eta*eta);
    dNdeta[8] = -0.5*(1-eps)*(2*eta);

    jacobianSurface(points, nbrnds, TotNbrNds, dNdeps, dNdeta);
}
void KERNEL::T3(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps,double eta)
{
    if (TotNbrNds != 3)
    {cout << "ERROR in KERNEL::T3. Total support points incompatible with element type" << endl; exit(0);}

    N.resize(TotNbrNds+1);
    vector<double> dNdeps(TotNbrNds+1);
    vector<double> dNdeta(TotNbrNds+1);
    
    N[1] = (1-eps-eta);
    N[2] = eps;
    N[3] = eta;

    dNdeps[1] = -1;
    dNdeta[1] = -1;
    
    dNdeps[2] = 1;
    dNdeta[2] = 0;
    
    dNdeps[3] = 0;
    dNdeta[3] = 1;

    jacobianSurface(points, nbrnds, TotNbrNds, dNdeps, dNdeta);
}
void  KERNEL::T6(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, double eps,double eta)
{
    if (TotNbrNds != 6)
    {cout << "ERROR in KERNEL::T6. Total support points incompatible with element type" << endl; exit(0);}

    N.resize(TotNbrNds+1);
    vector<double> dNdeps(TotNbrNds+1);
    vector<double> dNdeta(TotNbrNds+1);

    N[1] = (2*(1-eps-eta)-1)*(1-eps-eta);
    N[2] = eps*(2*eps - 1);
    N[3] = eta*(2*eta - 1);
    N[4] = 4*eps*(1-eps-eta);
    N[5] = 4*eps*eta;
    N[6] = 4*eta*(1-eps-eta);

    dNdeps[1] = (-(2*(1-eps-eta)-1))-2*(1-eps-eta);
    dNdeta[1] = (-(2*(1-eps-eta)-1))-2*(1-eps-eta);
    
    dNdeps[2] = 4*eps - 1;
    dNdeta[2] = 0;
    
    dNdeps[3] = 0;
    dNdeta[3] = 4*eta - 1;
    
    dNdeps[4] = -4*eps + 4*(1-eps-eta);
    dNdeta[4] = -4*eps;
    
    dNdeps[5] = 4*eta;
    dNdeta[5] = 4*eps;
    
    dNdeps[6] = -4*eta;
    dNdeta[6] = -4*eta + 4*(1-eps-eta);

    jacobianSurface(points, nbrnds, TotNbrNds, dNdeps, dNdeta);

}
void KERNEL::jacobianVolume(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, vector<double> dNdeps, vector<double> dNdeta, vector<double> dNdzta)
{
    vector<vector<double> > inverse_3X3(vector<vector<double> > A);
    double det_3X3(vector<vector<double> > mat);

    dNdx.resize(TotNbrNds+1);
    dNdy.resize(TotNbrNds+1);
    dNdz.resize(TotNbrNds+1);

    pnt.x = 0;
    pnt.y = 0;
    pnt.z = 0;

    double dxdeps = 0;
    double dxdeta = 0;
    double dxdzta = 0;
    double dydeps = 0;
    double dydeta = 0;
    double dydzta = 0;
    double dzdeps = 0;
    double dzdeta = 0;
    double dzdzta = 0;
    for (int i = 1; i <= TotNbrNds; i++)
    {
        pnt.x = pnt.x + N[i]* points[nbrnds[i]].x;
        pnt.y = pnt.y + N[i]* points[nbrnds[i]].y;
        pnt.z = pnt.z + N[i]* points[nbrnds[i]].z;
               
        dxdeps = dxdeps + dNdeps[i]*points[nbrnds[i]].x;
        dydeps = dydeps + dNdeps[i]*points[nbrnds[i]].y;
        dzdeps = dzdeps + dNdeps[i]*points[nbrnds[i]].z;

        dxdeta = dxdeta + dNdeta[i]*points[nbrnds[i]].x;
        dydeta = dydeta + dNdeta[i]*points[nbrnds[i]].y;
        dzdeta = dzdeta + dNdeta[i]*points[nbrnds[i]].z;

        dxdzta = dxdzta + dNdzta[i]*points[nbrnds[i]].x;
        dydzta = dydzta + dNdzta[i]*points[nbrnds[i]].y;
        dzdzta = dzdzta + dNdzta[i]*points[nbrnds[i]].z;
    }
    
    vector<vector<double> > jac(4,vector<double>(4));
    jac[1][1] = dxdeps;
    jac[1][2] = dydeps;
    jac[1][3] = dzdeps;

    jac[2][1] = dxdeta;
    jac[2][2] = dydeta;
    jac[2][3] = dzdeta;

    jac[3][1] = dxdzta;
    jac[3][2] = dydzta;
    jac[3][3] = dzdzta;

    det_jac = det_3X3(jac);
    vector<vector<double> > jacobian_inv(4, vector<double> (4));
    jacobian_inv = inverse_3X3(jac);

    for (int i = 1; i <= TotNbrNds; i++)
    {
        dNdx[i] = dNdeps[i]*jacobian_inv[1][1] + dNdeta[i]*jacobian_inv[1][2] + dNdzta[i]*jacobian_inv[1][3];
        dNdy[i] = dNdeps[i]*jacobian_inv[2][1] + dNdeta[i]*jacobian_inv[2][2] + dNdzta[i]*jacobian_inv[2][3];
        dNdz[i] = dNdeps[i]*jacobian_inv[3][1] + dNdeta[i]*jacobian_inv[3][2] + dNdzta[i]*jacobian_inv[3][3];
    }
    double sum = 0;
    for (int i = 1; i <= TotNbrNds; i++)
    {sum =sum + N[i];}
    if (sum <= 0.9)
    {cout << "ERROR in constructed Lagrange interpolants. Partition of Unity is Not satisfied. sum = " << sum << endl;exit(0);}
}
void KERNEL::jacobianSurface(vector<POINT> points, vector<int> nbrnds, int TotNbrNds, vector<double> dNdeps, vector<double> dNdeta)
{
    double det_2X2(vector<vector<double> > mat);

    dNdx.resize(TotNbrNds+1);
    dNdy.resize(TotNbrNds+1);
    dNdz.resize(TotNbrNds+1);

    pnt.x = 0;
    pnt.y = 0;
    pnt.z = 0;
    
    double dxdeps = 0;
    double dxdeta = 0;
    double dydeps = 0;
    double dydeta = 0;
    for (int i = 1; i <= TotNbrNds; i++)
    {
        pnt.x = pnt.x + N[i]* points[nbrnds[i]].x;
        pnt.y = pnt.y + N[i]* points[nbrnds[i]].y;
        pnt.z = pnt.z + N[i]* points[nbrnds[i]].z;

        dxdeps = dxdeps + dNdeps[i]*points[nbrnds[i]].x;
        dydeps = dydeps + dNdeps[i]*points[nbrnds[i]].y;

        dxdeta = dxdeta + dNdeta[i]*points[nbrnds[i]].x;
        dydeta = dydeta + dNdeta[i]*points[nbrnds[i]].y;
    }

    vector<vector<double> > jac(3,vector<double>(3));
    jac[1][1] = dxdeps;
    jac[2][1] = dydeps;

    jac[1][2] = dxdeta;
    jac[2][2] = dydeta;

    det_jac = det_2X2(jac);
    vector<vector<double> > jacobian_inv(3, vector<double> (3));
    jacobian_inv = inverse_2X2(jac);
    for (int i = 1; i <= TotNbrNds; i++)
    {
        dNdx[i] = dNdeps[i]*jacobian_inv[1][1] + dNdeta[i]*jacobian_inv[2][1];
        dNdy[i] = dNdeps[i]*jacobian_inv[1][2] + dNdeta[i]*jacobian_inv[2][2];
        dNdz[i] = 0;
    }
    double sum = 0;
    for (int i = 1; i <= TotNbrNds; i++)
    {sum =sum + N[i];}
    if (sum <= 0.9)
    {cout << "ERROR in constructed Lagrange interpolants. Partition of Unity is Not satisfied. sum = " << sum << endl;exit(0);}
}
KERNEL::~KERNEL()
{
    //destructor
}
