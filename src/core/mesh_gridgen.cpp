#include "../femwrks.h"
//----------------------------------------------------------------------
// FUNCTION DEFINITIONS FOR POINT CLASS
//----------------------------------------------------------------------
POINT::POINT()
{
}
POINT::POINT(double X , double Y , double Z)
{
    x = X;
    y = Y;
    z = Z;
}
POINT::POINT(double X , double Y)
{
    x = X;
    y = Y;
    z = 0;
}
POINT::POINT(double X)
{
    x = X;
    y = 0;
    z = 0;
}
double POINT::distance(POINT point)
{
    double dx = x-point.x;
    double dy = y-point.y;
    double dz = z-point.z;
    double distance = pow((dx*dx + dy*dy + dz*dz),0.5);
    return distance;
}
POINT::~POINT()
{
    //destructor
}
//----------------------------------------------------------------------
// FUNCTION DEFINITIONS FOR GRID DISCRITIZATION CLASS
//----------------------------------------------------------------------
RESOLUTION::RESOLUTION()
{
}
RESOLUTION::RESOLUTION(int nx, int ny, int nz)
{
	if (nx < 1 || ny < 1 || nz < 1)
	{
		cout << "ERROR: resolution can not be less than 1 point per axis" << endl;
		exit(0);
	}
	else
	{
        x = nx;
        y = ny;
        z = nz;
	}
}
RESOLUTION::RESOLUTION(int nx, int ny)
{
	if (nx < 1 || ny < 1)
	{
		cout << "ERROR: resolution can not be less than 1 point per axis" << endl;
		exit(0);
	}
	else
	{
		x = nx;
		y = ny;
		z = 1;
	}
}
RESOLUTION::RESOLUTION(int nx)
{
	if (nx < 1)
	{
		cout << "ERROR: resolution can not be less than 1 point per axis" << endl;
		exit(0);
	}
	else
	{
		x = nx;
		y = 1;
		z = 1;
	}
}
RESOLUTION::~RESOLUTION()
{
    //destructor
}

//----------------------------------------------------------------------
// FUNCTION DEFINITIONS FOR DOMAIN DOMAIN_SIZE CLASS
//----------------------------------------------------------------------
DOMAIN::DOMAIN()
{

}
DOMAIN::DOMAIN(double lx, double ly, double lz)
{
	if (lx < 0 || ly < 0 || lz < 0)
	{
		cout << "ERROR: domain size can not be less than 0" << endl;
		exit(0);
	}
	else
	{
		x = lx;
		y = ly;
        z = lz;
    }
}
DOMAIN::DOMAIN(double lx, double ly)
{
	if (lx < 0 || ly < 0)
	{
		cout << "ERROR: domain size can not be less than 0" << endl;
		exit(0);
	}
	else
	{
		x = lx;
		y = ly;
        z = 0;
    }
}
DOMAIN::DOMAIN(double lx)
{
	if (lx < 0)
	{
		cout << "ERROR: domain size can not be less than 0" << endl;
		exit(0);
	}
	else
	{
		x = lx;
		y = 0;
        z = 0;
    }
}
DOMAIN::~DOMAIN()
{
    //destructor
}

//----------------------------------------------------------------------
// FUNCTION DEFINITIONS FOR MESH CLASS
//----------------------------------------------------------------------
MESH::MESH() //constructor
{
    TotalLines= 0;
    TotalPoints = 0;
    TotalSurfaces = 0;
    TotalVolumes = 0;
    dim = 0;
}
MESH::MESH(POINT p, DOMAIN l, RESOLUTION divss, string elemtype)
{
    cout << "Generating grid..." << endl;
    if (divss.z == 1)
    {
        if (divss.y > 1)
        {
            dim = 2;
            generateGridNodes2D(p, l, divss);
            if (elemtype == "T3")
            {generateGridTri(divss);}
            else if (elemtype == "Q4")
            {generateGridQuad(divss);}
            else
            {cout << "ERROR: unrecognized 2D element type. Allowable values are T3 or Q4" << endl; exit(0);}
        }
        else if(divss.y == 1)
        {
            dim = 1;
            generateGridNodes1D(p, l, divss);
        }
    }
    else if (divss.z > 1)
    {
        dim = 3;
        generateGridNodes3D(p,l,divss);
        if (elemtype == "T4")
        {generateGridTet(p, l, divss);}
        else if (elemtype == "H8")
        {generateGridHex(p, l, divss);}
        else
        {cout << "ERROR: unrecognized 3D element type. Allowable values are T4 or H8" << endl; exit(0);}
    }
    cout << "Total Nodes = " << TotalPoints << endl;
    cout << "Total Lines = " << TotalLines << endl;
    cout << "Total Surfaces = " << TotalSurfaces << endl;
    cout << "Total Volumes = " << TotalVolumes << endl;
    cout << "Mesh read successfuly..." << endl;
}

//------------------------------------------------------------------
//FUNCTIONS RELATED TO GENERATING GRID POINTS
//------------------------------------------------------------------
void MESH::generateGridNodes1D(POINT p, DOMAIN s, RESOLUTION divs)
{
    points.resize(divs.x+1);
    double spacing_x = s.x/double(divs.x-1);
    points[1].x = p.x;
    points[1].y = 0;
    points[1].z = 0;
    int l = 1;
    for (int i = 1; i <= divs.x; i++)
    {
        points[l+1].x = points[l].x + spacing_x;
        points[l+1].y = 0;
        points[l+1].z = 0;
        l = l+1;
    }
    TotalPoints = l-1;
}
void MESH::generateGridNodes2D(POINT p1, DOMAIN s, RESOLUTION divs)
{
    points.resize((divs.x+1)*(divs.y+1));
    double spacing_x = s.x/double(divs.x-1);
    double spacing_y = s.y/double(divs.y-1);
    points[1].x = p1.x;
    points[1].y = p1.y;
    points[1].z = 0;
    int l = 1;
    for (int j = 1; j <= divs.y; j++)
    {
        for (int i = 1; i <= divs.x; i++)
        {
            points[l+1].x = points[l].x + spacing_x;
            points[l+1].y = points[l].y;
            points[l+1].z = 0;
            l = l+1;
        }
        points[l].x = points[1].x;
        points[l].y = points[l].y + spacing_y;
        points[l].z = 0;
    }
    TotalPoints = l-1;
}
void MESH::generateGridNodes3D(POINT p1, DOMAIN s, RESOLUTION divs)
{
    points.resize((divs.x+1)*(divs.y+1)*(divs.z+1));
    double spacing_x = s.x/double(divs.x-1);
    double spacing_y = s.y/double(divs.y-1);
    double spacing_z = s.z/double(divs.z-1);
    points[1].x = p1.x;
    points[1].y = p1.y;
    points[1].z = p1.z;
    int l = 1;
    for (int k = 1; k <= divs.z; k++)
    {
        for (int j = 1; j <= divs.y; j++)
        {
            for (int i = 1; i <= divs.x; i++)
            {
                points[l+1].x = points[l].x + spacing_x;
                points[l+1].y = points[l].y;
                points[l+1].z = points[l].z;
                l++;
            }
            points[l].x = points[1].x;
            points[l].y = points[l].y + spacing_y;
            points[l].z = points[l].z;
        }
        points[l].x = p1.x;
        points[l].y = p1.y;
        points[l].z = points[l].z+spacing_z;
    }
    TotalPoints = l-1;
    cout << "Total constructed grid points = " << TotalPoints << endl;
}
//------------------------------------------------------------------
//FUNCTIONS RELATED TO GENERATING GRID ELEMENTS
//------------------------------------------------------------------
void MESH::generateGridQuad(RESOLUTION divs)
{
    int elems_x = divs.x-1;
    int elems_y = divs.y-1;
    dim = 2;
    ElementOrder = 1;
    SurfaceShape = 2;
    SurfaceNds = 4;
    
    int elem_num = elems_x * elems_y;
    surfaces.resize(elem_num+1);
    int node_num = (elems_x+1)*(elems_y+1);
    int node_n[node_num+1];
    for (int i = 1; i <= node_num; i++)
    {node_n[i] = i;}

    int k = 1;
    int l = 1;
    for (int j = 1; j<= elems_y; j++)
    {
        for (int i = 1; i <= elems_x; i++)
        {
            surfaces[l].resize(SurfaceNds+1);
            surfaces[l][1] = node_n[k];
            surfaces[l][2] = node_n[k+1];
            surfaces[l][3] = node_n[k+elems_x+2];
            surfaces[l][4] = node_n[k+elems_x+1];
            k = k+1;
            l = l+1;
        }
        k = k+1;
    }
    TotalLines = 0;
    TotalSurfaces = l-1;
    TotalVolumes = 0;
}
void MESH::generateGridTri(RESOLUTION divs)
{
    int elems_x = divs.x-1;
    int elems_y = divs.y-1;
    dim = 2;
    ElementOrder = 1;
    SurfaceShape = 1;
    SurfaceNds = 3;
                
    int elem_num = 2*elems_x*elems_y;
    surfaces.resize(elem_num+1);
    int node_num = (elems_x+1)*(elems_y+1);
    int node_n[node_num+1];
    for (int i = 1; i <= node_num; i++)
    {node_n[i] = i;}

    int k = 1;
    int l = 1;
    for (int j = 1; j<= elems_y; j++)
    {
        for (int i = 1; i <= elems_x; i++)
        {
            if ((mod(i,2) == 0 && mod(j,2) != 0) || (mod(i,2) != 0 && mod(j,2) == 0))
            {

                surfaces[l].resize(SurfaceNds+1);
                surfaces[l][1] = node_n[k];
                surfaces[l][2] = node_n[k+1];
                surfaces[l][3] = node_n[k+elems_x+1];

                surfaces[l+1].resize(SurfaceNds+1);
                surfaces[l+1][1] = node_n[k+1];
                surfaces[l+1][2] = node_n[k+elems_x+2];
                surfaces[l+1][3] = node_n[k+elems_x+1];
            }
            else if ((mod(i,2) != 0 && mod(j,2) != 0) || (mod(i,2) == 0 && mod(j,2)== 0))
            {
                surfaces[l].resize(SurfaceNds+1);
                surfaces[l][1] = node_n[k];
                surfaces[l][2] = node_n[k+1];
                surfaces[l][3] = node_n[k+elems_x+2];

                surfaces[l+1].resize(SurfaceNds+1);
                surfaces[l+1][1] = node_n[k];
                surfaces[l+1][2] = node_n[k+elems_x+2];
                surfaces[l+1][3] = node_n[k+elems_x+1];
            }
            k = k+1;
            l = l+2;
        }
        k = k+1;
    }
    TotalLines = 0;
    TotalSurfaces = l-1;
    TotalVolumes = 0;
}
void MESH::generateGridTri2(RESOLUTION divs)
{
    int elems_x = divs.x-1;
    int elems_y = divs.y-1;
    dim = 2;
    ElementOrder = 1;
    SurfaceShape = 1;
    SurfaceNds = 3;
    
    int elem_num = 2*elems_x * elems_y;
    surfaces.resize(elem_num+1);
    int node_num = (elems_x+1)*(elems_y+1);
    int node_n[node_num+1];
    for (int i = 1; i <= node_num; i++)
    {node_n[i] = i;}

    int k = 1;
    int l = 1;
    for (int j = 1; j<= elems_y; j++)
    {
        for (int i = 1; i <= elems_x; i++)
        {
            surfaces[l].resize(SurfaceNds+1);
            surfaces[l][1] = node_n[k];
            surfaces[l][2] = node_n[k+1];
            surfaces[l][3] = node_n[k+elems_x+1];

            surfaces[l+1].resize(SurfaceNds+1);
            surfaces[l+1][1] = node_n[k+1];
            surfaces[l+1][2] = node_n[k+elems_x+2];
            surfaces[l+1][3] = node_n[k+elems_x+1];
            l = l+2;
            k = k+1;
        }
        k = k+1;
    }
    TotalLines = 0;
    TotalSurfaces = elem_num;
    TotalVolumes = 0;
}
void MESH::generateGridTet(POINT p, DOMAIN size, RESOLUTION divs)
{
    dim = 3;
    ElementOrder = 1;
    VolumeShape = 1;
    SurfaceShape = 1;
    VolumeNds = 4;
    SurfaceNds = 3;
    
    vector<vector<vector<int> > > TotalPoints(divs.x+1, vector<vector<int> >(divs.y+1, vector<int>(divs.z+1)));
    TotalLines = 0;
    TotalVolumes = 5*(divs.x-1)*(divs.y-1)*(divs.z-1);
    volumes.resize(TotalVolumes+1);
    int e = 1;

    int l = 1;
    for (int k = 1; k <= divs.z; k++)
    {
        for (int j = 1; j <= divs.y; j++)
        {
            for (int i = 1; i <= divs.x; i++)
            {
                TotalPoints[i][j][k] = l;
                l++;
            }
        }
    }

    for (int k = 1; k <= divs.z-1; k++)
    {
        for (int j = 1; j <= divs.y-1; j++)
        {
            for (int i = 1; i <= divs.x-1; i++)
            {
                volumes[e].resize(VolumeNds+1);
                volumes[e][1] = TotalPoints[i][j][k];
                volumes[e][2] = TotalPoints[i+1][j][k];
                volumes[e][3] = TotalPoints[i][j+1][k];
                volumes[e][4] = TotalPoints[i][j][k+1];

                volumes[e+1].resize(VolumeNds+1);
                volumes[e+1][1] = TotalPoints[i+1][j][k];
                volumes[e+1][2] = TotalPoints[i+1][j+1][k];
                volumes[e+1][3] = TotalPoints[i][j+1][k];
                volumes[e+1][4] = TotalPoints[i+1][j+1][k+1];

                volumes[e+2].resize(VolumeNds+1);
                volumes[e+2][1] = TotalPoints[i][j][k+1];
                volumes[e+2][2] = TotalPoints[i][j+1][k];
                volumes[e+2][3] = TotalPoints[i][j+1][k+1];
                volumes[e+2][4] = TotalPoints[i+1][j+1][k+1];

                volumes[e+3].resize(VolumeNds+1);
                volumes[e+3][1] = TotalPoints[i][j][k+1];
                volumes[e+3][2] = TotalPoints[i+1][j][k];
                volumes[e+3][3] = TotalPoints[i+1][j+1][k+1];
                volumes[e+3][4] = TotalPoints[i+1][j][k+1];

                volumes[e+4].resize(VolumeNds+1);
                volumes[e+4][1] = TotalPoints[i][j+1][k];
                volumes[e+4][2] = TotalPoints[i+1][j][k];
                volumes[e+4][3] = TotalPoints[i][j][k+1];
                volumes[e+4][4] = TotalPoints[i+1][j+1][k+1];

                e = e+5;
            }
        }
    }
    TotalSurfaces = TotalVolumes*4;
    surfaces.resize(TotalSurfaces+1);
    int s = 1;
    for (int i = 1; i <= TotalVolumes; i++)
    {
        surfaces[s].resize(SurfaceNds+1);
        surfaces[s][1] = volumes[i][1];
        surfaces[s][2] = volumes[i][2];
        surfaces[s][3] = volumes[i][3];

        surfaces[s+1].resize(SurfaceNds+1);
        surfaces[s+1][1] = volumes[i][1];
        surfaces[s+1][2] = volumes[i][2];
        surfaces[s+1][3] = volumes[i][4];

        surfaces[s+2].resize(SurfaceNds+1);
        surfaces[s+2][1] = volumes[i][2];
        surfaces[s+2][2] = volumes[i][3];
        surfaces[s+2][3] = volumes[i][4];

        surfaces[s+3].resize(SurfaceNds+1);
        surfaces[s+3][1] = volumes[i][1];
        surfaces[s+3][2] = volumes[i][3];
        surfaces[s+3][3] = volumes[i][4];
        s = s + 4;
    }
}
void MESH::generateGridHex(POINT p, DOMAIN size, RESOLUTION divs)
{
    dim = 3;
    ElementOrder = 1;
    VolumeShape = 2;
    SurfaceShape = 2;
    VolumeNds = 8;
    SurfaceNds = 4;
    
    vector<vector<vector<int> > > TotalPoints(divs.x+1, vector<vector<int> >(divs.y+1, vector<int>(divs.z+1)));
    TotalLines = 0;
    TotalVolumes = (divs.x-1)*(divs.y-1)*(divs.z-1);
    volumes.resize(TotalVolumes+1);
    int e = 1;

    int l = 1;
    for (int k = 1; k <= divs.z; k++)
    {
        for (int j = 1; j <= divs.y; j++)
        {
            for (int i = 1; i <= divs.x; i++)
            {
                TotalPoints[i][j][k] = l;
                l++;
            }
        }
    }

    for (int k = 1; k <= divs.z-1; k++)
    {
        for (int j = 1; j <= divs.y-1; j++)
        {
            for (int i = 1; i <= divs.x-1; i++)
            {
                volumes[e].resize(VolumeNds+1);
                volumes[e][1] = TotalPoints[i][j][k];
                volumes[e][2] = TotalPoints[i+1][j][k];
                volumes[e][3] = TotalPoints[i+1][j+1][k];
                volumes[e][4] = TotalPoints[i][j+1][k];
                volumes[e][5] = TotalPoints[i][j][k+1];
                volumes[e][6] = TotalPoints[i+1][j][k+1];
                volumes[e][7] = TotalPoints[i+1][j+1][k+1];
                volumes[e][8] = TotalPoints[i][j+1][k+1];
                e = e+1;
            }
        }
    }
    TotalSurfaces = TotalVolumes*6;
    surfaces.resize(TotalSurfaces+1);
    int s = 1;
    for (int i = 1; i <= TotalVolumes; i++)
    {
        surfaces[s].resize(SurfaceNds+1);
        surfaces[s][1] = volumes[i][1];
        surfaces[s][2] = volumes[i][2];
        surfaces[s][3] = volumes[i][3];
		surfaces[s][4] = volumes[i][4];

        surfaces[s+1].resize(SurfaceNds+1);
        surfaces[s+1][1] = volumes[i][5];
        surfaces[s+1][2] = volumes[i][6];
        surfaces[s+1][3] = volumes[i][7];
        surfaces[s+1][4] = volumes[i][8];

        surfaces[s+2].resize(SurfaceNds+1);
        surfaces[s+2][1] = volumes[i][2];
        surfaces[s+2][2] = volumes[i][6];
        surfaces[s+2][3] = volumes[i][7];
		surfaces[s+2][4] = volumes[i][3];

        surfaces[s+3].resize(SurfaceNds+1);
        surfaces[s+3][1] = volumes[i][1];
        surfaces[s+3][2] = volumes[i][5];
        surfaces[s+3][3] = volumes[i][8];
        surfaces[s+3][4] = volumes[i][4];

		surfaces[s+4].resize(SurfaceNds+1);
        surfaces[s+4][1] = volumes[i][1];
        surfaces[s+4][2] = volumes[i][2];
        surfaces[s+4][3] = volumes[i][6];
        surfaces[s+4][4] = volumes[i][5];

        surfaces[s+5].resize(SurfaceNds+1);
        surfaces[s+5][1] = volumes[i][4];
        surfaces[s+5][2] = volumes[i][3];
        surfaces[s+5][3] = volumes[i][7];
        surfaces[s+5][4] = volumes[i][8];
        
        s = s + 6;
    }
}
MESH::~MESH()
{
    //destructor
}

