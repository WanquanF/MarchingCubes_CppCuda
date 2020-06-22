# marching_cubes_cuda
cmake_minimum_required (VERSION 2.8)
project(Test)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_CXX_FLAGS "-std=c++11 -w -D_REENTRANT -fPIC -pipe -g -fopenmp ${CMAKE_CXX_FLAGS}")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -w")

include_directories(${CMAKE_SOURCE_DIR})

#CUDA
find_package(CUDA REQUIRED)
set(CUDA_ARCH "30 35 50 52 60 61" CACHE STRING "CUDA architectures to compile for")
string(REPLACE " " ";" CUDA_ARCH_LIST ${CUDA_ARCH})
foreach(a ${CUDA_ARCH_LIST})
  set(CUDA_ARCH_OPTIONS "${CUDA_ARCH_OPTIONS} -gencode arch=compute_${a},code=compute_${a}")
endforeach()
include_directories(${CUDA_INCLUDE_DIRS})
#include_directories("/usr/local/cuda-8.0/targets/x86_64-linux/include")
#link_directories("/usr/local/cuda-8.0/targets/x86_64-linux/lib")
include_directories("/usr/local/cuda-8.0/include")
link_directories("/usr/local/cuda-8.0/lib64")

#Gl, GLOG
link_directories("/usr/lib/x86_64-linux-gnu")
#dlib, boost_system, protobuf, blas, lapack
link_directories("/usr/lib")


include_directories("/home/wanquan/eigen3")
include_directories("/disk3/wanquan/meshsr/include/libigl/include")


cuda_add_executable(marching_cubes_cuda
    cuda_inc.h
    cuda_func.h
    common.h
    mc_glue.h
    mc_glue.cpp
    mc_glue.cu
    mc.cpp
)

target_link_libraries(marching_cubes_cuda
        -lpthread
        -lcuda 
        -lcudart 
        -lcublas 
        -lcusolver 
)




