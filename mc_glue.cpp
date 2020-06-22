#include "mc_glue.h"


mc_glue::mc_glue(int x_grid_num,
            int y_grid_num,
            int z_grid_num,
            std::string sdf_path_in,
            std::string mesh_path_out)
{

    
    all_points_record.clear();
    all_faces.clear();
	fTargetValue=0.0;
    
    sdf_file_ = sdf_path_in;
    mesh_path_out_ = mesh_path_out;
    NX=x_grid_num;
    NY=y_grid_num;
    NZ=z_grid_num;
	
	
	number_record_=(int*)malloc(2*sizeof(int));
	number_record_[0]=0;
	number_record_[1]=0;
	
	cuda_init();
	set_length_and_read_voxel();
	
    
}


void mc_glue::set_length_and_read_voxel()
{
    sdf_readin_=(float*)malloc(NX*NY*NZ*sizeof(float));
    
    FILE *fid;
    fid = fopen(sdf_file_.c_str(),"rb");
    for(int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ;k++)
            {
                float oc;
                fread(&oc,sizeof(float),1,fid);
                sdf_readin_[i*NY*NZ+j*NZ+k]=oc;
            }
        }
    }
    fclose(fid);
    
    edge_point_state=(int*)malloc(NX*NY*NZ*3*sizeof(int));
    for (int i=0;i<NX;i++)
    {
        for (int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ;k++)
            {
                for(int l=0;l<3;l++)
                {
                    edge_point_state[i*NY*NZ*3+j*NZ*3+k*3+l]=-1;
                }
            }
        }
    }
    
    cuda_get_sdf_values();
}


float mc_glue::fGetOffset(float fValue1, float fValue2, float fValueDesired)
{
    double fDelta = fValue2 - fValue1;

    if(fDelta == 0.0)
    {
        return 0.5;
    }
    return (fValueDesired - fValue1)/fDelta;
}


void mc_glue::vMarchCube1(float fX, float fY, float fZ)
{
    //cout<<"dealing with  "<<fX<<" "<<fY<<" "<<fZ<<endl;
    int iCorner, iVertex, iVertexTest, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
    float fOffset;
    float afCubeValue[8];
    float asEdgeVertex[12][3];
    //coor_vector asEdgeNorm[12];

    for(iVertex = 0; iVertex < 8; iVertex++)
    {
        afCubeValue[iVertex] =(float)sdf_readin_[(int)(fX + a2fVertexOffset[iVertex][0])*NY*NZ+(int)(fY + a2fVertexOffset[iVertex][1])*NZ+(int)(fZ + a2fVertexOffset[iVertex][2])];
        //cout<<a2fVertexOffset[iVertex][0]<<","<<a2fVertexOffset[iVertex][1]<<","<<a2fVertexOffset[iVertex][2]<<":"<<afCubeValue[iVertex]<<endl;
    }

    iFlagIndex = 0;
    for(iVertexTest = 0; iVertexTest < 8; iVertexTest++)
    {
        if(afCubeValue[iVertexTest] < fTargetValue)
            iFlagIndex |= 1<<iVertexTest;
    }

    iEdgeFlags = aiCubeEdgeFlags[iFlagIndex];

    if(iEdgeFlags == 0)
    {
        return;
    }
    for(iEdge = 0; iEdge < 12; iEdge++)
    {
        if(iEdgeFlags & (1<<iEdge))
        {
            fOffset = fGetOffset(afCubeValue[ a2iEdgeConnection[iEdge][0] ],
                    afCubeValue[ a2iEdgeConnection[iEdge][1] ], fTargetValue);


            asEdgeVertex[iEdge][0] = fX + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][0]  +  fOffset * a2fEdgeDirection[iEdge][0]) ;//* fScale;
            asEdgeVertex[iEdge][1] = fY + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][1]  +  fOffset * a2fEdgeDirection[iEdge][1]) ;//* fScale;
            asEdgeVertex[iEdge][2] = fZ + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][2]  +  fOffset * a2fEdgeDirection[iEdge][2]) ;//* fScale;

            //vGetNormal(asEdgeNorm[iEdge], asEdgeVertex[iEdge].fX, asEdgeVertex[iEdge].fY, asEdgeVertex[iEdge].fZ);
        }
    }


    for(iTriangle = 0; iTriangle < 5; iTriangle++)
    {
        //std::cout<<"There is a new triangle."<<std::endl;
        if(a2iTriangleConnectionTable[iFlagIndex][3*iTriangle] < 0)
            break;

        Eigen::Vector3i af;

        for(iCorner = 0; iCorner < 3; iCorner++)
        {
            iVertex = a2iTriangleConnectionTable[iFlagIndex][3*iTriangle+iCorner];
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
            if(edge_point_state[basex*NY*NZ*3+basey*NZ*3+basez*3+direction]==-1)
            {
                edge_point_state[basex*NY*NZ*3+basey*NZ*3+basez*3+direction]=all_points_record.size();
                Eigen::Vector3f av;
                av<<
                     asEdgeVertex[iVertex][0],asEdgeVertex[iVertex][1],asEdgeVertex[iVertex][2];
                all_points_record.push_back(av);
                af(iCorner)=all_points_record.size()-1;
            }else
            {
                af(iCorner)=edge_point_state[basex*NY*NZ*3+basey*NZ*3+basez*3+direction];
            }
        }
        all_faces.push_back(af);
    }
}

void mc_glue::mc_get_mesh()
{
    clock_t s,e;
    s=clock();
    
    
    for(int i=0;i<NX-1;i++)
    {
        for(int j=0;j<NY-1;j++)
        {
            for(int k=0;k<NZ-1;k++)
            {
                vMarchCube1(i,j,k);
            }
        }
    }

    
    std::ofstream ouf;
    ouf.open(mesh_path_out_);
    for (int v =0; v<all_points_record.size();v++)
    {
        ouf<<"v "<<all_points_record[v](0)<<" "<<all_points_record[v](1)<<" "<<all_points_record[v](2)<<std::endl;
    }
    for (int f =0; f<all_faces.size();f++)
    {
        ouf<<"f "<<all_faces[f](2)+1<<" "<<all_faces[f](1)+1<<" "<<all_faces[f](0)+1<<std::endl;
    }
    ouf.close();
    std::cout<<"The result mesh has been saved."<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;
}

void mc_glue::mc_get_mesh_cuda()
{
    mc_get_mesh_on_gpu();
    
    cudaMemcpy(number_record_,d_number_record_,sizeof(int)*2,cudaMemcpyDeviceToHost);
    
//    std::cout<<number_record_[0]<<" "<<number_record_[1]<<std::endl;
    convert_ijkd_to_pindex_on_gpu();
    
    cudaMemcpy(number_record_,d_number_record_,sizeof(int)*2,cudaMemcpyDeviceToHost);
    faces_index_ = (int*)malloc(number_record_[1]*3*sizeof(int));
    points_coor_ = (float*)malloc(number_record_[0]*3*sizeof(float));;
    
    cudaMemcpy(faces_index_,d_faces_index_,sizeof(int)*number_record_[1]*3,cudaMemcpyDeviceToHost);
    cudaMemcpy(points_coor_,d_points_coor_,sizeof(float)*number_record_[0]*3,cudaMemcpyDeviceToHost);
    
    std::ofstream ouf;
    ouf.open(mesh_path_out_);
    for (int v =0; v<number_record_[0];v++)
    {
        ouf<<"v "<<points_coor_[v*3]<<" "<<points_coor_[v*3+1]<<" "<<points_coor_[v*3+2]<<std::endl;
    }
    for (int f =0; f<number_record_[1];f++)
    {
        ouf<<"f "<<faces_index_[f*3]+1<<" "<<faces_index_[f*3+1]+1<<" "<<faces_index_[f*3+2]+1<<std::endl;
    }
    ouf.close();
    std::cout<<"The result mesh has been saved."<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;
    
}















