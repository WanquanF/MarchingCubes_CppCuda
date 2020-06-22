#include "mc_glue.h"
#include "cuda_func.h"

void mc_glue::cuda_init()
{


    cudaFree(d_a2fVertexOffset_);
    CHECK(cudaMalloc((float**)&d_a2fVertexOffset_,sizeof(float)*8*3));
    cudaMemcpy(d_a2fVertexOffset_,a2fVertexOffset,sizeof(float)*8*3,cudaMemcpyHostToDevice);
    
    cudaFree(d_a2iEdgeConnection_);
    CHECK(cudaMalloc((int**)&d_a2iEdgeConnection_,sizeof(int)*12*2));
    cudaMemcpy(d_a2iEdgeConnection_,a2iEdgeConnection,sizeof(int)*12*2,cudaMemcpyHostToDevice);
    
    cudaFree(d_a2fEdgeDirection_);
    CHECK(cudaMalloc((float**)&d_a2fEdgeDirection_,sizeof(float)*12*3));
    cudaMemcpy(d_a2fEdgeDirection_,a2fEdgeDirection,sizeof(float)*12*3,cudaMemcpyHostToDevice);
    
    cudaFree(d_aiCubeEdgeFlags_);
    CHECK(cudaMalloc((int**)&d_aiCubeEdgeFlags_,sizeof(int)*256));
    cudaMemcpy(d_aiCubeEdgeFlags_,aiCubeEdgeFlags,sizeof(int)*256,cudaMemcpyHostToDevice);
    
    cudaFree(d_a2iTriangleConnectionTable_);
    CHECK(cudaMalloc((int**)&d_a2iTriangleConnectionTable_,sizeof(int)*256*16));
    cudaMemcpy(d_a2iTriangleConnectionTable_,a2iTriangleConnectionTable,sizeof(int)*256*16,cudaMemcpyHostToDevice);
    
    cudaFree(d_number_record_);
    CHECK(cudaMalloc((int**)&d_number_record_,sizeof(int)*2));
    cudaMemcpy(d_number_record_,number_record_,sizeof(int)*2,cudaMemcpyHostToDevice);
    
}


void mc_glue::cuda_get_sdf_values()
{
    cudaFree(d_sdf_);
    CHECK(cudaMalloc((float**)&d_sdf_,sizeof(float)*NX*NY*NZ));
    cudaMemcpy(d_sdf_,sdf_readin_,sizeof(float)*NX*NY*NZ,cudaMemcpyHostToDevice);
    
    cudaFree(d_edge_point_state_);
    CHECK(cudaMalloc((int**)&d_edge_point_state_,sizeof(int)*NX*NY*NZ*3));
    cudaMemcpy(d_edge_point_state_,edge_point_state,sizeof(int)*NX*NY*NZ*3,cudaMemcpyHostToDevice);
    
    
    cudaFree(d_points_coor_);
    CHECK(cudaMalloc((float**)&d_points_coor_,sizeof(float)*3*(int)(NX*NY*NZ*12*0.05)));
    
    cudaFree(d_faces_index_);
    CHECK(cudaMalloc((int**)&d_faces_index_,sizeof(int)*3*(int)(NX*NY*NZ*5*0.05)));
    
    cudaFree(d_faces_ijkd_);
    CHECK(cudaMalloc((int**)&d_faces_ijkd_,sizeof(int)*3*4*(int)(NX*NY*NZ*5*0.05)));
}


void mc_glue::mc_get_mesh_on_gpu()
{
    int n_thread=NX*NY*NZ;
    d_mc_get_mesh_on_gpu<<<GET_BLOCKS(n_thread),CUDA_NUM_THREADS>>>(n_thread,d_sdf_,d_edge_point_state_,NX,NY,NZ, fTargetValue, d_number_record_,d_a2fVertexOffset_,d_a2iEdgeConnection_,d_a2fEdgeDirection_,d_aiCubeEdgeFlags_,d_a2iTriangleConnectionTable_,d_points_coor_,d_faces_ijkd_,d_faces_index_);
}


void mc_glue::convert_ijkd_to_pindex_on_gpu()
{
    int n_thread=number_record_[1];
    d_conver_ijkd_to_pindex<<<GET_BLOCKS(n_thread),CUDA_NUM_THREADS>>>(n_thread,NX,NY,NZ, d_edge_point_state_, d_faces_ijkd_, d_faces_index_, d_number_record_);
}

