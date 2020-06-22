import os
import sys

os.system('./build/marching_cubes_cuda 205 205 205 ./sample/cell_rot_0_low.sdf ./sample/cell_rot_0_low_rec_with_cuda.obj 1')

os.system('./build/marching_cubes_cuda 205 205 205 ./sample/cell_rot_0_low.sdf ./sample/cell_rot_0_low_rec_without_cuda.obj 0')
