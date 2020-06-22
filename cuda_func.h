#ifndef CUDA_FUNC_H
#define CUDA_FUNC_H

#include "cuda_inc.h"


inline __host__ __device__ float d_fGetOffset(float fValue1, float fValue2, float fValueDesired)
{
    double fDelta = fValue2 - fValue1;

    if(fDelta == 0.0)
    {
        return 0.5;
    }
    return (fValueDesired - fValue1)/fDelta;
}


__global__ void d_mc_get_mesh_on_gpu(int nthread, float* d_sdf, int *d_edge_point_state, int NX, int NY, int NZ, float fTargetValue, int *d_number_record, float* d_a2fVertexOffset, int *d_a2iEdgeConnection, float *d_a2fEdgeDirection, int* d_aiCubeEdgeFlags, int *d_a2iTriangleConnectionTable, float *d_points_coor, int *d_faces_ijkd, int *d_faces_index)
{
    CUDA_KERNEL_LOOP(index, nthread)
    {
        int i,j,k;
        i=index/(NY*NZ);
        j=(index-i*NY*NZ)/NZ;
        k=index-i*NY*NZ-j*NZ;
        if( i<NX-1 && j<NY-1 && k<NZ-1)
//        if(true)
        {
            float fX, fY, fZ;
            fX = i;
            fY = j;
            fZ = k;
            
            int iCorner, iVertex, iVertexTest, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
            float fOffset;
            float afCubeValue[8];
            float asEdgeVertex[12][3];
            
            for(iVertex = 0; iVertex < 8; iVertex++)
            {
                int index_to_use=(int)(fX + d_a2fVertexOffset[iVertex*3+0])*NY*NZ+(int)(fY + d_a2fVertexOffset[iVertex*3+1])*NZ+(int)(fZ + d_a2fVertexOffset[iVertex*3+2]);
                afCubeValue[iVertex] =d_sdf[index_to_use];
            }
            
            iFlagIndex = 0;
            for(iVertexTest = 0; iVertexTest < 8; iVertexTest++)
            {
                if(afCubeValue[iVertexTest] < fTargetValue)
                {
                    iFlagIndex |= 1<<iVertexTest;
                }
            }

            iEdgeFlags = d_aiCubeEdgeFlags[iFlagIndex];
            
            
            
            if(iEdgeFlags != 0)
            {
                
                for(iEdge = 0; iEdge < 12; iEdge++)
                {
                    if(iEdgeFlags & (1<<iEdge))
                    {
                        fOffset = d_fGetOffset(afCubeValue[ d_a2iEdgeConnection[iEdge*2+0] ], afCubeValue[ d_a2iEdgeConnection[iEdge*2+1] ], fTargetValue);
                        asEdgeVertex[iEdge][0] = fX + (d_a2fVertexOffset[d_a2iEdgeConnection[iEdge*2+0]*3+0] + fOffset * d_a2fEdgeDirection[iEdge*3+0]);
                        asEdgeVertex[iEdge][1] = fY + (d_a2fVertexOffset[d_a2iEdgeConnection[iEdge*2+0]*3+1] + fOffset * d_a2fEdgeDirection[iEdge*3+1]);
                        asEdgeVertex[iEdge][2] = fZ + (d_a2fVertexOffset[d_a2iEdgeConnection[iEdge*2+0]*3+2] + fOffset * d_a2fEdgeDirection[iEdge*3+2]);
                    }
                }
            
            
            
                bool is_ivertex_new[12];
                for (int iin =0;iin<12;iin++)
                {
                    is_ivertex_new[iin]=true;
                }
                for(iTriangle = 0; iTriangle < 5; iTriangle++)
                {
                    
                    if(d_a2iTriangleConnectionTable[iFlagIndex*16+3*iTriangle] < 0)
                    {
                        break;
                    }
                
                    int face_id = atomicAdd(&d_number_record[1],1);

                    for(iCorner = 0; iCorner < 3; iCorner++)
                    {
                        iVertex = d_a2iTriangleConnectionTable[iFlagIndex*16+3*iTriangle+iCorner];
                        int basex,basey,basez;
                        int direction;
                        if(iVertex==0)
                        {
                            basex=fX;
                            basey=fY;
                            basez=fZ;
                            direction=0;
                        }else if(iVertex==1)
                        {
                            basex=fX+1;
                            basey=fY;
                            basez=fZ;
                            direction=1;
                        }else if(iVertex==2)
                        {
                            basex=fX;
                            basey=fY+1;
                            basez=fZ;
                            direction=0;
                        }else if(iVertex==3)
                        {
                            basex=fX;
                            basey=fY;
                            basez=fZ;
                            direction=1;
                        }else if(iVertex==4)
                        {
                            basex=fX;
                            basey=fY;
                            basez=fZ+1;
                            direction=0;
                        }else if(iVertex==5)
                        {
                            basex=fX+1;
                            basey=fY;
                            basez=fZ+1;
                            direction=1;
                        }else if(iVertex==6)
                        {
                            basex=fX;
                            basey=fY+1;
                            basez=fZ+1;
                            direction=0;
                        }else if(iVertex==7)
                        {
                            basex=fX;
                            basey=fY;
                            basez=fZ+1;
                            direction=1;
                        }else if(iVertex==8)
                        {
                            basex=fX;
                            basey=fY;
                            basez=fZ;
                            direction=2;
                        }else if(iVertex==9)
                        {
                            basex=fX+1;
                            basey=fY;
                            basez=fZ;
                            direction=2;
                        }else if(iVertex==10)
                        {
                            basex=fX+1;
                            basey=fY+1;
                            basez=fZ;
                            direction=2;
                        }else if(iVertex==11)
                        {
                            basex=fX;
                            basey=fY+1;
                            basez=fZ;
                            direction=2;
                        }
                        
                        int vert_id=-1;
                        if (is_ivertex_new[iVertex]==false)
                        {
//                            vert_id = atomicMax(&d_edge_point_state[basex*NY*NZ*3+basey*NZ*3+basez*3+direction],-10);
                        }
                        else if (iVertex==0 || iVertex==3 || iVertex==8)
                        {
                            vert_id = atomicAdd(&d_number_record[0],1);
                            d_points_coor[vert_id*3+0]=asEdgeVertex[iVertex][0];
                            d_points_coor[vert_id*3+1]=asEdgeVertex[iVertex][1];
                            d_points_coor[vert_id*3+2]=asEdgeVertex[iVertex][2];
                            atomicExch(&d_edge_point_state[basex*NY*NZ*3+basey*NZ*3+basez*3+direction],vert_id);
                            is_ivertex_new[iVertex]=false;
                        }
                        d_faces_ijkd[face_id*3*4+iCorner*4+0]=basex;
                        d_faces_ijkd[face_id*3*4+iCorner*4+1]=basey;
                        d_faces_ijkd[face_id*3*4+iCorner*4+2]=basez;
                        d_faces_ijkd[face_id*3*4+iCorner*4+3]=direction;
                    }
                }
            }
        }
    }
}


__global__ void d_conver_ijkd_to_pindex(int nthread, int NX, int NY, int NZ, int *d_edge_point_state, int *d_faces_ijkd, int *d_faces_index, int *d_number_record)
{
    CUDA_KERNEL_LOOP(index, nthread)
    {
        for (int pid = 0; pid<3; pid++)
        {
            int i = d_faces_ijkd[index*3*4+pid*4+0];
            int j = d_faces_ijkd[index*3*4+pid*4+1];
            int k = d_faces_ijkd[index*3*4+pid*4+2];
            int d = d_faces_ijkd[index*3*4+pid*4+3];
            d_faces_index[index*3+(2-pid)]=d_edge_point_state[i*NY*NZ*3+j*NZ*3+k*3+d];
        }
    }
}




#endif // CUDA_FUNC_H
































