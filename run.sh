#!/bin/bash
echo "*****************************"
echo "3D FINITE ELEMENT SOLVER"
echo "*****************************"
echo "Generating makefile..."
mkdir out
mkdir obj
mkdir pics
rm -r out/*
rm -r pics/*
rm obj/main.o
rm run
rm a.out
rm makefile.mak
g++ src/makefile.cpp
./a.out
echo "Compiling source code..."
make -f makefile.mak
echo "Running executable..."
./run
echo "Cleaning..."
rm run
rm makefile.mak
rm a.out
rm obj/main.o
echo "Script complete..."
