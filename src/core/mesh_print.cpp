#include "../femwrks.h"

static string output_dir = "out";

//-------------------------------------------------------------------------------
// FUNCTIONS FOR PRINTING MESH IN VTK FORMAT
//-------------------------------------------------------------------------------
void MESH::printMeshVTK(vector<double> U, string filename, int t)
{
    cout << "Printing mesh data to " << filename << "..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_vtk/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    stringstream tt;
    tt << t;

    if (t>=1 && t <= 9)
    {fname =  direc + filename+ "00"+tt.str()+".vtk";}
    else if (t >= 10 && t <= 99)
    {fname =  direc + filename+ "0"+tt.str()+".vtk";}
    else
    {fname =  direc + filename+ tt.str()+".vtk";}
    int sum = 0;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        int num=0;
        if (SurfaceShape == 1)
        {num = 4;}
        else if (SurfaceShape == 2)
        {num = 5;}
        else
        {cout << "ERROR2: unrecognized element type" << endl;}
        sum = sum + num;
    }
    ofstream outfile1(fname.c_str());
    outfile1 << "# vtk DataFile Version 1.0" << endl;
    outfile1 << "Cube example" << endl;
    outfile1 << "ASCII" << endl;
    outfile1 << endl;
    outfile1 <<"DATASET POLYDATA" << endl;
    outfile1 << "POINTS " << TotalPoints << " double" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}
    outfile1 << endl;
    outfile1 << "POLYGONS " << TotalSurfaces << "\t" << sum << endl;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        outfile1 << SurfaceNds/ElementOrder << "\t";
        for (int j = 1; j <= SurfaceNds/ElementOrder; j++)
        {outfile1 << (surfaces[i][j]-1)<< "\t";}
        outfile1 << endl;
    }
    outfile1 << endl;
    outfile1 <<"POINT_DATA " << TotalPoints << endl;
    outfile1 << "SCALARS myscalars double"<< endl;
    outfile1 << "LOOKUP_TABLE custom_table" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << U[i] << endl;}
    outfile1.close();
}
void MESH::printMeshVTK(vector<int> U, string filename, int t)
{
    cout << "Printing mesh data to " << filename << "..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_vtk/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    stringstream tt;
    tt << t;

    if (t>=1 && t <= 9)
    {fname =  direc + filename+ "00"+tt.str()+".vtk";}
    else if (t >= 10 && t <= 99)
    {fname =  direc + filename+ "0"+tt.str()+".vtk";}
    else
    {fname =  direc + filename+ tt.str()+".vtk";}
    int sum = 0;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        int num=0;
        if (SurfaceShape == 1)
        {num = 4;}
        else if (SurfaceShape == 2)
        {num = 5;}
        else
        {cout << "ERROR2: unrecognized element type" << endl;}
        sum = sum + num;
    }
    ofstream outfile1(fname.c_str());
    outfile1 << "# vtk DataFile Version 1.0" << endl;
    outfile1 << "Cube example" << endl;
    outfile1 << "ASCII" << endl;
    outfile1 << endl;
    outfile1 <<"DATASET POLYDATA" << endl;
    outfile1 << "POINTS " << TotalPoints << " double" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}
    outfile1 << endl;
    outfile1 << "POLYGONS " << TotalSurfaces << "\t" << sum << endl;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        outfile1 << SurfaceNds/ElementOrder << "\t";
        for (int j = 1; j <= SurfaceNds/ElementOrder; j++)
        {outfile1 << (surfaces[i][j]-1)<< "\t";}
        outfile1 << endl;
    }
    outfile1 << endl;
    outfile1 <<"POINT_DATA " << TotalPoints << endl;
    outfile1 << "SCALARS myscalars double"<< endl;
    outfile1 << "LOOKUP_TABLE custom_table" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << U[i] << endl;}
    outfile1.close();
}

void MESH::printVectorsVTK(vector<double> Vx, vector<double> Vy, vector<double> Vz, string filename, int t)
{
    cout << "Printing mesh data to " << filename << "..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_vtk/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    stringstream tt;
    tt << t;

    if (t>=1 && t <= 9)
    {fname = direc + filename+ "00"+tt.str()+".vtk";}
    else if (t >= 10 && t <= 99)
    {fname = direc + filename+ "0"+tt.str()+".vtk";}
    else
    {fname = direc + filename+ tt.str()+".vtk";}

    vector<double> V(TotalPoints+1);
    for (int i = 1; i <= TotalPoints; i++)
    {V[i] = pow((pow(Vx[i],2)+pow(Vy[i],2)+pow(Vz[i],2)),0.5);}

    int sum = 0;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        int num=0;
        if (SurfaceShape == 1)
        {num = 4;}
        else if (SurfaceShape == 2)
        {num = 5;}
        else
        {cout << "ERROR2: unrecognized element type. SurfaceShape = " << SurfaceShape << endl;}
        sum = sum + num;
    }
    ofstream outfile1(fname.c_str());
    outfile1 << "# vtk DataFile Version 1.0" << endl;
    outfile1 << "Cube example" << endl;
    outfile1 << "ASCII" << endl;
    outfile1 << endl;
    outfile1 <<"DATASET POLYDATA" << endl;
    outfile1 << "POINTS " << TotalPoints << " double" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}
    outfile1 << endl;
    outfile1 << "POLYGONS " << TotalSurfaces << "\t" << sum << endl;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        outfile1 << SurfaceNds/ElementOrder << "\t";
        for (int j = 1; j <= SurfaceNds/ElementOrder; j++)
        {outfile1 << (surfaces[i][j]-1)<< "\t";}
        outfile1 << endl;
    }
    outfile1 << endl;
    outfile1 <<"POINT_DATA " << TotalPoints << endl;
    outfile1 << "SCALARS myscalars double"<< endl;
    outfile1 << "LOOKUP_TABLE custom_table" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << V[i] << endl;}
    outfile1 << endl;
    outfile1 << "VECTORS vectors double" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << Vx[i] << "\t" << Vy[i] <<"\t" << Vz[i] << endl;}
    outfile1.close();
}
void MESH::printTXT(vector<double> U, string filename, int t)
{
    cout << "Printing solution to file..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_txt/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    stringstream tt;
    tt << t;
    if (t>=1 && t <= 9)
    {fname = direc + filename+ "00"+tt.str()+".txt";}
    else if (t >= 10 && t <= 99)
    {fname = direc + filename+ "0"+tt.str()+".txt";}
    else
    {fname = direc + filename+ tt.str()+".txt";}
    ofstream outfile(fname.c_str());
    for (int i = 1; i <= TotalPoints; i++)
    {outfile << points[i].x << "\t" << points[i].y<<  "\t"  << points[i].z << "\t" << U[i] << endl;}
    outfile.close();
}
void printElements(vector<vector<int> > surfaces, int SurfaceNds, int TotalSurfaces, string filename)
{
    cout << "Printing mesh data to " << filename << "..." << endl;
    string direc = output_dir+"/gmtry/";
    mkdir(direc.c_str(), 0777);
    string fname;
    fname = direc + filename + ".txt";
    ofstream outfile1(fname.c_str());
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        for (int j = 1; j <= SurfaceNds; j++)
        {outfile1 << (surfaces[i][j]-1)<< "\t";}
        outfile1 << endl;
    }
    outfile1.close();
}
void pltctrl(vector<POINT> points, int TotalPoints, int n_t, int print_frequency)
{
    string direc = output_dir+"/gmtry/";
    mkdir(direc.c_str(), 0777);

    string fname;
    double xmin = points[1].x;
    double xmax = points[1].x;
    double ymin = points[1].y;
    double ymax = points[1].y;
    for (int i = 1; i <= TotalPoints; i++)
    {
        if (points[i].x <= xmin)
        {xmin =  points[i].x;}
        if (points[i].x >= xmax)
        {xmax =  points[i].x;}
        if (points[i].y <= ymin)
        {ymin =  points[i].y;}
        if (points[i].y >= ymax)
        {ymax =  points[i].y;}
    }
    fname = output_dir+"/gmtry/"+"pltctrl.txt";
    ofstream outfile(fname.c_str());
    outfile << n_t << "\t" << print_frequency << "\t" << xmin << "\t" <<  xmax << "\t" << ymin << "\t" << ymax << endl;
    outfile.close();
}

void MESH::print()
{
    string direc;
    direc = output_dir+"/" + "meshoutput/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    fname = direc + "nodes.txt";
    ofstream outfile1(fname.c_str());
    outfile1 << TotalPoints << "\t" << 3 << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}
        
        
    fname = direc + "surfaces.txt";
    ofstream outfile2(fname.c_str());
    outfile2 << TotalSurfaces << "\t" << SurfaceNds << "\t" << SurfaceShape << "\t" << ElementOrder << endl;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
		for (int j = 1; j <= SurfaceNds; j++)
		{outfile2 << surfaces[i][j] << "\t";}
		outfile2 << endl;
	}   
	
	    fname = direc + "volumes.txt";
    ofstream outfile3(fname.c_str());
    outfile3 << TotalVolumes << "\t" << VolumeNds << "\t" << VolumeShape << "\t" << ElementOrder << endl;
    for (int i = 1; i <= TotalVolumes; i++)
    {
		for (int j = 1; j <= VolumeNds; j++)
		{outfile3 << volumes[i][j] << "\t";}
		outfile3 << endl;
	}   
	
    outfile1.close();
    outfile2.close();
    outfile3.close();
}
