#!/bin/bash

cd Butane
octave -q run_butane.m
cd ..

cd Cflow
octave -q run_Cflow.m
cd ..

cd Diatomic
octave -q run_diatomic.m
cd ..

cd KA
octave -q run_KA.m
cd ..

cd LJ
octave -q run_LJ.m
cd ..

cd Toluene
octave -q run_toluene.m
cd ..

