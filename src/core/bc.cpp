#include "../femwrks.h"

BOUNDARY_CONDITION::BOUNDARY_CONDITION() //constructor
{
    total = 0;
    points.resize(1);
	values.resize(1);
}
void BOUNDARY_CONDITION::assignDBC(MESH mesh, string salomeface_fname, double value)
{
    cout << "Assigning surface direchlet boundary conditions..." << endl;
    if (salomeface_fname == "")
    {cout << "ERROR: direchlet boundary conditions file name is invalid..." << endl;exit(0);}
    else
    {
        MESH boundary(salomeface_fname);
        if (boundary.TotalPoints == 0)
        {cout << "ERROR: direchlet boundary conditions points are empty..." << endl;exit(0);}
        else
        {
            vector<int> pnt_num(boundary.TotalPoints+1);
            for (int i = 1; i <= boundary.TotalPoints; i++)
            {
                for (int j = 1; j <= mesh.TotalPoints; j++)
                {
                    double d = boundary.points[i].distance(mesh.points[j]);
                    if (d <= 1e-5)
                    {pnt_num[i] = j;}
                }
            }
            for (int i = 1; i <= boundary.TotalPoints; i++)
            {
				points.push_back(1);
				values.push_back(1);
				total++;
				points[total] = pnt_num[i];
				values[total] = value;

			}
            cout << "SUCCESS: direchlet boundary conditions read successfully..." << endl;
        }
    }
}
void BOUNDARY_CONDITION::assignDBC(MESH mesh, string edge, double location, double value)
{
    cout << "Assigning grid direchlet boundary conditions..." << endl;
    vector<int> ndnum(1);
    int m = 0;
    double tol = 0.01;
    if (edge == "x")
    {
        for (int i = 1; i <= mesh.TotalPoints; i++)
        {
            if (mesh.points[i].x >= location-tol*abs(location) && mesh.points[i].x <= location+tol*abs(location))
            {
                m++;
                ndnum.push_back(1);
                ndnum[m] = i;
            }
        }
    }
    if (edge == "y")
    {
        for (int i = 1; i <= mesh.TotalPoints; i++)
        {
            if (mesh.points[i].y >= location-tol*abs(location) && mesh.points[i].y <= location+tol*abs(location))
            {
                m++;
                ndnum.push_back(1);
                ndnum[m] = i;
            }
        }
    }
    if (edge == "z")
    {
        for (int i = 1; i <= mesh.TotalPoints; i++)
        {
            if (mesh.points[i].z >= location-tol*abs(location) && mesh.points[i].z <= location+tol*abs(location))
            {
                m++;
                ndnum.push_back(1);
                ndnum[m] = i;
            }
        }
    }

    if (m == 0)
    {cout << "ERROR: direchlet boundary conditions nodes are empty..." << endl;exit(0);}
    else
    {
        for (int i = 1; i <= m; i++)
        {
			total++;
			points.push_back(1);
			values.push_back(1);
			points[total] = ndnum[i];
			values[total] = value;
		}
        cout << "SUCCESS: direchlet boundary conditions assigned successfully..." << endl;
        cout << "Total direchlet boundary nodes = " << m << endl;
    }
}

BOUNDARY_CONDITION::~BOUNDARY_CONDITION()
{
    //destructor
}

vector<double> setVector(vector<double> U, BOUNDARY_CONDITION direchletBC)
{
    for (int r = 1; r <= direchletBC.total; r++)
    {U[direchletBC.points[r]] = direchletBC.values[r];}
    return U;
}

vector<double> setVector(int TotalPoints, double val)
{
	vector<double> vec(TotalPoints+1);
	for (int i = 1; i <= TotalPoints; i++)
	{vec[i] = val;}
	return vec;
}
