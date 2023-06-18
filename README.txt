-----------------------------------------------------------------------------------
FemWrks
----------------------------------------------------------------------------------
FemWrks is a 2D/3D finite element solver written in C++ and is used to solve
the general transport equation, incompressible navier-stokes equations and can
be used for stress/deformation analysis. It can use first and second order
elements and can solve the finite element equations using element-by-element
technique without the need to form global matrices.
---------------------------------------------------------------------------------
FemWrks has been tested to run on Debian and FreeBSD systems. To run FemWrks perform
the following steps:

1. On Debian, ensure the required packages are installed:
sudo apt install bash gcc g++ make cmake libomp-dev openmpi-bin libgomp1 build-essential python3 python3-numpy python3-scipy python3-matplotlib

On FreeBSD, ensure the following packages are installed:
sudo pkg install gcc openmpi openmpi3 gmake cmake python python3 py38-numpy py38-scipy py38-matplotlib

2. Inside the FemWrks folder, edit permissions to shell scripts:
chmod u+x run.sh cln.sh genvid.sh

3. Copy a case study from the examples folder and paste inside the main FemWrks folder. Rename it to "main.cpp"

4. Run compiler script:
./run.sh

or on FreeBSD:
bash run.sh

5. All the simulation results are generated and stored inside the "out" folder. 
To view the results data in "VTK" format, use any software that open vtk files such as ParaView, VisIt, Salome-platform

6. Data generated in "TXT" can be used to generate jpg images of the simulation using the included python script:
python3 pltmsh.py <scale prmtr> <show msh flag> <var 1> <var 2>

where the first parameter controls the scaling of the plot, the second is 1/0 which shows/hides the mesh, the third parameter is the variable you want plotted.

example:
python3 pltmsh.py 5 1 U

7. To export to video from the generated jpg images do the following:

./genvid.sh U

where U is the name of the variable that was used to generate the simulation pics in step 5


