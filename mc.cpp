#include <iostream>
using namespace std;
#include <ctime>
#include <time.h>

#include "mc_glue.h"

int main(int argc,char **argv)
{
    int x_grid_num=atoi(argv[1]);
    cout<<"x_grid_num = "<<x_grid_num<<endl;
    int y_grid_num=atoi(argv[2]);
    cout<<"y_grid_num = "<<y_grid_num<<endl;
    int z_grid_num = atoi(argv[3]);
    cout<<"z_grid_num = "<<z_grid_num<<endl;
    string sdf_path_in = argv[4];
    cout<<"sdf_path_in = "<<sdf_path_in<<endl;
    string mesh_path_out = argv[5];
    cout<<"mesh_path_out = "<<mesh_path_out<<endl;
    int with_cuda = atoi(argv[6]);
    if (with_cuda >0)
    {
        cout<<"with_cuda = Yes "<<endl;
    }else
    {
        cout<<"with_cuda = No "<<endl;
    }
    
    mc_glue* mger=new mc_glue(x_grid_num,
                              y_grid_num,
                              z_grid_num, sdf_path_in, mesh_path_out);
    
    if(with_cuda>0)
    {
        mger->mc_get_mesh_cuda();   
    }else
    {
        mger->mc_get_mesh();
    }
    
    delete mger;
    mger=NULL; 
    
    
    
    return 0;
}
