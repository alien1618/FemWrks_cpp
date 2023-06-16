#include<iostream>
#include<fstream>

using namespace std;

//TO INSTALL IN DEBIAN DO THIS:
//TO SUCCESSFULLY COMPILE YOU NEED GCC AND OPENMP INSTALLED
// sudo apt-get install gcc libgomp1

int main()
{
    ofstream file("makefile.mak");
    file << "CPP = c++" << endl;
    file << "CC = gcc" << endl;
    file << "RES = " << endl;
    file << "VPATH = src: src/core/: src/slvrs: " << endl;
    file << "OBJDIR = obj" << endl;

    //--------------------------------------------------------------------------
    //CORE FILES
    //--------------------------------------------------------------------------
    file << "OBJ =  $(OBJDIR)/mesh_gridgen.o ";
    file << "$(OBJDIR)/mesh_read.o ";
    file << "$(OBJDIR)/mesh_print.o ";
    file << "$(OBJDIR)/quadrature.o ";
    file << "$(OBJDIR)/eqslvrs.o ";
    file << "$(OBJDIR)/matrixops.o ";
    file << "$(OBJDIR)/globalmatrices.o ";
    file << "$(OBJDIR)/bc.o ";
    file << "$(OBJDIR)/interpolation.o ";
    file << "$(OBJDIR)/matrixfree.o ";

    //--------------------------------------------------------------------------
    //SOLVERS
    //--------------------------------------------------------------------------
    file << "$(OBJDIR)/slv_transport.o ";
    file << "$(OBJDIR)/slv_navierstokes.o ";
    file << "$(OBJDIR)/slv_elasticity.o ";

    //--------------------------------------------------------------------------
    //MAIN
    //--------------------------------------------------------------------------
    file << "$(OBJDIR)/main.o " << " $(RES)" << endl;

    file << "LINKOBJ = $(OBJ)" << endl;
    file << "BIN  = run" << endl;
    file << "CXXFLAGS = $(CXXINCS)" << endl;
    file <<"CFLAGS = $(INCS)" << endl;

    file <<"RM = rm -f" << endl;

    file <<".PHONY: all clean clean-custom" << endl;

    file <<"all: run" << endl;

    file <<"clean: clean-custom" << endl;
    file <<"	${RM} $(OBJ) $(BIN)" << endl;

    file <<"$(BIN): $(OBJ)" << endl;
    file <<"	$(CPP) $(LINKOBJ) -o \"run\" $(LIBS) -std=c++17 -fopenmp -Wall -O3" << endl;

    file << "$(OBJDIR)/%.o: %.cpp" << endl;
    file <<"	$(CPP) -c $(CXXFLAGS) $< -o $@ -std=c++17 -fopenmp -Wall -O3" << endl;
    
    file.close();
    return 0;
}
