#!/bin/bash

module purge
module load gcc/11.1.0
module load cmake/3.20.2
module load mpich/3.3.0

source /home/aurelien.vadrot/TOOLS/miniconda3/etc/profile.d/conda.sh
conda activate env_smarties
export CC=`which gcc`
export CXX=`which g++`


cd build
rm -rf *
cmake ..
make -j
#cd ../apps/LESGO_M
#make clean
#make
