#include"nMDS_R.h"
using namespace std;

extern "C"{
void nMDS_R(int *number_of_points,int *profile_length,double *data,int *embedding_dimensions,int *number_of_iterations, double *output_positions,int *random_seed,int *distance_measure)
{
  srand(*random_seed);
  nMDS nMDS_Instance(*number_of_points,*profile_length,data,*embedding_dimensions,*number_of_iterations,*distance_measure);
  nMDS_Instance.Output(output_positions);
}
}
