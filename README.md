# MarchingCubes_CppCuda
The C++ (with cuda) version of the classic Marching Cubes algorithm.

# Description
This repo contains the code for the classic Marching Cubes algorithm in C++&Cuda. The input is a SDF tensor with shape (dim-x, dim-y, dim-z). The output is the result reconstructed mesh by the MC algorithm.  

## Prerequisite Installation
    1. eigen3
    2. cuda-8.0
    3. Python3
    
## How to use the code: 
Clone the repo:
``` bash
    git clone https://github.com/WanquanF/MarchingCubes_CppCuda.git
```
Build the code with the given script:
``` bash
    python run.py
```
Run the testing sample:
``` bash
    python test.py
```
## Other details:
The code has two kinds of implementations of the MC algorithm. One implementation uses the cpu all the time while the other one uses the gpu for parallel. After you build the code by running the "run.py" script, the exe file is in the "./build" folder, named "marching_cubes_cuda". The command line for using it is like this:
``` bash
    ./build/marching_cubes_cuda  dim-x  dim-y  dim-z  input-sdf-path  output-mesh-path  if-using-cuda
```
The "dim-x", "dim-y" and "dim-z" are int numbers showing the shape of the input sdf, and the "if-using-cuda" is a 0/1 (0 means running without cuda, 1 means running with cuda) value deciding which version of implementation you use. The "input-sdf-path" is a bin file of the sdf value tensor. The "output-mesh-path" should be in the obj format.
