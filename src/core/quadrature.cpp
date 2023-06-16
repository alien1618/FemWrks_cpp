#include "../femwrks.h"

void GAUSS::generatePnts(int dim, int elem_type, int level)
{
    //-------------------------------------------------
    // generate gauss points per element
    //-------------------------------------------------
    switch(dim)
    {
        case 3:
        {
            switch(level)
            {
                case 1:
                {
                    switch(elem_type)
                    {
                        case 1: //tetrahedrals
                        {
                            tetra(1);
                            element_integpoints = 1;
                        }break;
                        case 2: //hexahedrons
                        {
                            hexa(2);
                            element_integpoints = int(pow(2.0,3.0));
                        }break;
                    }
                }break;
                case 2:
                {
                    switch (elem_type)
                    {
                        case 1: //tetrahedrals
                        {
                            tetra(4);
                            element_integpoints = 4;
                        }break;
                        case 2: //hexahedrons
                        {
                            hexa(2);
                            element_integpoints = int(pow(2,3.0));
                        }break;
                    }
                }break;
            }
        }break;
        case 2:
        {
            switch (level)
            {
                case 1:
                {
                    switch(elem_type)
                    {
                        case 1: //triangular
                        {
                            tri(1);
                            element_integpoints = 1;
                        }break;
                        case 2: //quads
                        {
                            quad(2);
                            element_integpoints = 4;
                        }break;
                    }
                }break;
                case 2:
                {
                    switch (elem_type)
                    {
                        case 1: //triangular
                        {
                            tri(3);
                            element_integpoints = 3;
                        }break;
                        case 2: //quads
                        {
                            quad(2);
                            element_integpoints = 4;
                        }break;
                    }
                }break;
                case 3:
                {
                    switch (elem_type)
                    {
                        case 1: //triangular
                        {
                            tri(7);
                            element_integpoints = 7;
                        }break;
                        case 2: //quads
                        {
                            quad(6);
                            element_integpoints = 36;
                        }break;
                    }
                }break;
            }
        }break;
        case 1:
        {
            element_integpoints = 6;
            line(element_integpoints);
        }break;
    }
}
void GAUSS::hexa(int num_pts)
{
    double a1, a2,a3, w1,w2,w3;
    vector<double> eps_pnt(num_pts+1);
    vector<double> w_pnt(num_pts+1);
    switch (num_pts)
    {
        case 1:
        {
            a1 = 0.000000000000000;
            w1 = 2.000000000000000;
            eps_pnt[1] = a1;
            w_pnt[1] = w1;
        }break;
        case 2:
        {
            a1 = 0.577350269189626;
            w1 = 1.000000000000000;
            eps_pnt[1] = a1;
            eps_pnt[2] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w1;
        }break;
        case 3:
        {
            a1 = 0.774596669241483;
            a2 = 0.000000000000000;
            w1 = 0.555555555555556;
            w2 = 0.888888888888889;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w1;
        }break;
        case 4:
        {
            a1 = 0.861136311594053;
            a2 = 0.339981043584856;
            w1 = 0.347854845137454;
            w2 = 0.652145154862546;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = -a2;
            eps_pnt[4] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w2;
            w_pnt[4] = w1;
        }break;
        case 5:
        {
            a1 = 0.906179845938664;
            a2 = 0.538469310105683;
            a3 = 0.000000000000000;
            w1 = 0.236926885056189;
            w2 = 0.478628670499366;
            w3 = 0.568888888888889;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = a3;
            eps_pnt[4] = -a2;
            eps_pnt[5] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w3;
            w_pnt[4] = w2;
            w_pnt[5] = w1;
        }break;
        case 6:
        {
            a1 = 0.932469514203152;
            a2 = 0.661209386466265;
            a3 = 0.238619186083197;
            w1 = 0.171324492379170;
            w2 = 0.360761573048139;
            w3 = 0.467913934572691;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = a3;
            eps_pnt[4] = -a3;
            eps_pnt[5] = -a2;
            eps_pnt[6] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w3;
            w_pnt[4] = w3;
            w_pnt[5] = w2;
            w_pnt[6] = w1;
        }break;
    }

    int integ_pts = num_pts*num_pts*num_pts;
    eps.resize(integ_pts+1);
    eta.resize(integ_pts+1);
    zta.resize(integ_pts+1);
    w.resize(integ_pts+2);
    int m = 0;
    for (int k = 1; k <= num_pts; k++)
    {
        for (int j = 1; j <= num_pts; j++)
        {
            for (int i = 1; i <= num_pts; i++)
            {
                m = m + 1;
                eps[m] = eps_pnt[i];
                eta[m] = eps_pnt[j];
                zta[m] = eps_pnt[k];
                w[m] = w_pnt[i]*w_pnt[j]*w_pnt[k];
            }
        }
    }
}
void GAUSS::tetra(int num_pts)
{
    eps.resize(num_pts+1);
    eta.resize(num_pts+1);
    zta.resize(num_pts+1);
    w.resize(num_pts+1);
    double a1, a2, w1;
    switch(num_pts)
    {
        case 1:
        {
            eps[1] = 0.250;
            eta[1] =  0.250;
            zta[1] =  0.250;
            w[1] = 1;
        }break;
        case 4:
        {
            a1 = 0.58541020;
            a2 = 0.13819660;
            w1 = 0.250;

            eps[1] = a1;
            eta[1] = a2;
            zta[1] = a2;
            w[1] = w1;

            eps[2] = a2;
            eta[2] = a1;
            zta[2] = a2;
            w[2] = w1;

            eps[3] = a2;
            eta[3] = a2;
            zta[3] = a1;
            w[3] = w1;

            eps[4] = a2;
            eta[4] = a2;
            zta[4] = a2;
            w[4] = w1;
        }break;
    }
}
void GAUSS::quad(int num_pts)
{
    double a1, a2,a3, w1,w2,w3;
    vector<double> eps_pnt(num_pts+1);
    vector<double> w_pnt(num_pts+1);
    switch(num_pts)
    {
        case 1:
        {
            a1 = 0.000000000000000;
            w1 = 2.000000000000000;
            eps_pnt[1] = a1;
            w_pnt[1] = w1;
        }break;
        case 2:
        {
            a1 = 0.577350269189626;
            w1 = 1.000000000000000;
            eps_pnt[1] = a1;
            eps_pnt[2] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w1;
        }break;
        case 3:
        {
            a1 = 0.774596669241483;
            a2 = 0.000000000000000;
            w1 = 0.555555555555556;
            w2 = 0.888888888888889;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w1;
        }break;
        case 4:
        {
            a1 = 0.861136311594053;
            a2 = 0.339981043584856;
            w1 = 0.347854845137454;
            w2 = 0.652145154862546;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = -a2;
            eps_pnt[4] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w2;
            w_pnt[4] = w1;
        }break;
        case 5:
        {
            a1 = 0.906179845938664;
            a2 = 0.538469310105683;
            a3 = 0.000000000000000;
            w1 = 0.236926885056189;
            w2 = 0.478628670499366;
            w3 = 0.568888888888889;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = a3;
            eps_pnt[4] = -a2;
            eps_pnt[5] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w3;
            w_pnt[4] = w2;
            w_pnt[5] = w1;
        }break;
        case 6:
        {
            a1 = 0.932469514203152;
            a2 = 0.661209386466265;
            a3 = 0.238619186083197;
            w1 = 0.171324492379170;
            w2 = 0.360761573048139;
            w3 = 0.467913934572691;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = a3;
            eps_pnt[4] = -a3;
            eps_pnt[5] = -a2;
            eps_pnt[6] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w3;
            w_pnt[4] = w3;
            w_pnt[5] = w2;
            w_pnt[6] = w1;
        }break;
    }

    int integ_pts = num_pts*num_pts;
    eps.resize(integ_pts+1);
    eta.resize(integ_pts+1);
    zta.resize(integ_pts+1);
    w.resize(integ_pts+1);
    int m = 0;
    for (int j = 1; j <= num_pts; j++)
    {
        for (int i = 1; i <= num_pts; i++)
        {
            m = m + 1;
            eps[m] = eps_pnt[i];
            eta[m] = eps_pnt[j];
            zta[m] = 1;
            w[m] = w_pnt[i]*w_pnt[j];
        }
    }

}
void GAUSS::tri(int num_pts)
{
    int integ_pts = num_pts*num_pts;
    eps.resize(integ_pts+1);
    eta.resize(integ_pts+1);
    zta.resize(integ_pts+1);
    w.resize(integ_pts+2);
    double a1, a2,a3, a4, a5, w1,w2,w3;
    switch(num_pts)
    {
        case 1:
        {
            eps[1] = 0.333333333333;
            eta[1] =  0.333333333333;
            zta[1] =  1;
            w[1] = 1;
        }break;
        case 3:
        {
            a1 = 0.1666666666667;
            a2 = 0.6666666666667;
            w1 = 0.1666666666667;
            eps[1] = a1;
            eps[2] = a2;
            eps[3] = a1;
            eta[1] = a1;
            eta[2] = a1;
            eta[3] = a2;
            zta[1] = 1;
            zta[2] = 1;
            zta[3] = 1;
            w[1] = w1;
            w[2] = w1;
            w[3] = w1;
        }break;
        case 7:
        {
            a1 = 0.1012865073235;
            a2 = 0.7974269853531;
            a3 = 0.4701420641051;
            a4 = 0.0597158717898;
            a5 = 0.3333333333333;
            w1 = 0.0629695902724;
            w2 = 0.0661970763942;
            w3 = 0.1125000000000;

            eps[1] = a1;
            eps[2] = a2;
            eps[3] = a1;
            eps[4] = a3;
            eps[5] = a3;
            eps[6] = a4;
            eps[7] = a5;
            eta[1] = a1;
            eta[2] = a1;
            eta[3] = a2;
            eta[4] = a4;
            eta[5] = a3;
            eta[6] = a3;
            eta[7] = a5;
            zta[1] = 1;
            zta[2] = 1;
            zta[3] = 1;
            zta[4] = 1;
            zta[5] = 1;
            zta[6] = 1;
            zta[7] = 1;
            w[1] = w1;
            w[2] = w1;
            w[3] = w1;
            w[4] = w2;
            w[5] = w2;
            w[6] = w2;
            w[7] = w3;
        } break;
    }
}
void GAUSS::line(int num_pts)
{
    double a1, a2,a3, w1,w2,w3;
    vector<double> eps_pnt(num_pts+1);
    vector<double> w_pnt(num_pts+1);
    switch(num_pts)
    {
        case 1:
        {
            a1 = 0.000000000000000;
            w1 = 2.000000000000000;
            eps_pnt[1] = a1;
            w_pnt[1] = w1;
        }break;
        case 2:
        {
            a1 = 0.577350269189626;
            w1 = 1.000000000000000;
            eps_pnt[1] = a1;
            eps_pnt[2] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w1;
        }break;
        case 3:
        {
            a1 = 0.774596669241483;
            a2 = 0.000000000000000;
            w1 = 0.555555555555556;
            w2 = 0.888888888888889;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w1;
        }break;
        case 4:
        {
            a1 = 0.861136311594053;
            a2 = 0.339981043584856;
            w1 = 0.347854845137454;
            w2 = 0.652145154862546;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = -a2;
            eps_pnt[4] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w2;
            w_pnt[4] = w1;
        }break;
        case 5:
        {
            a1 = 0.906179845938664;
            a2 = 0.538469310105683;
            a3 = 0.000000000000000;
            w1 = 0.236926885056189;
            w2 = 0.478628670499366;
            w3 = 0.568888888888889;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = a3;
            eps_pnt[4] = -a2;
            eps_pnt[5] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w3;
            w_pnt[4] = w2;
            w_pnt[5] = w1;
        }break;
        case 6:
        {
            a1 = 0.932469514203152;
            a2 = 0.661209386466265;
            a3 = 0.238619186083197;
            w1 = 0.171324492379170;
            w2 = 0.360761573048139;
            w3 = 0.467913934572691;
            eps_pnt[1] = a1;
            eps_pnt[2] = a2;
            eps_pnt[3] = a3;
            eps_pnt[4] = -a3;
            eps_pnt[5] = -a2;
            eps_pnt[6] = -a1;
            w_pnt[1] = w1;
            w_pnt[2] = w2;
            w_pnt[3] = w3;
            w_pnt[4] = w3;
            w_pnt[5] = w2;
            w_pnt[6] = w1;
        }break;
    }
    int integ_pts = num_pts;
    eps.resize(integ_pts+1);
    eta.resize(integ_pts+1);
    w.resize(integ_pts+1);
    for (int i = 1; i <= num_pts; i++)
    {
        eps[i] = eps_pnt[i];
        w[i] = w_pnt[i];
    }
}
GAUSS::~GAUSS()
{
    //destructor
}

