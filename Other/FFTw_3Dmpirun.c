#include "header.h"

int main(int argc, char **argv){

  MPI_Init(&argc,&argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &process_number);

  const int   lenPositionArray = 330;      // Number of cells per dimension
  const float CellSize         = 10.0;     // Size of the cell in units h^-1 Mpc 
  const int   TotalCellNumber  = lenPositionArray*lenPositionArray*lenPositionArray;  // Total number of cells in the box. 
  const       ptrdiff_t n0     = lenPositionArray, n1 = lenPositionArray, n2 = lenPositionArray;

  struct floatOutputArray{
    float floatArray[lenPositionArray*lenPositionArray*lenPositionArray];
  }; 

  ptrdiff_t alloc_local, local_n0, local_0_start, i, j, k;

  NyquistWavenumber  = pi/CellSize;                   // k = 2*pi x frequency        

  kInterval          = (float) 2*(pi)*pow(lenPositionArray, -1)*pow(CellSize, -1);

  int   NyquistIndex = lenPositionArray/2 + 1;

  fftw_mpi_init();

  /* get local data size and allocate */
  alloc_local        = fftw_mpi_local_size_3d(n0, n1, n2, MPI_COMM_WORLD, &local_n0, &local_0_start);

  float temp[alloc_local];

  in                 = fftw_alloc_complex(alloc_local);

  if(process_rank==0){
    printf("\n Master rank:  %d", process_rank);

    Array            = malloc((8*n0*n1*n2)*sizeof(*Array));
    root_dir         = malloc(200*sizeof(*root_dir));
    filepath         = malloc(200*sizeof(*filepath));

    sprintf(root_dir,"/disk1/mjw/ZADE_MockRun/");
    sprintf(filepath, "%sNGP_zadeMocks_w1_ByCellNumber_delta_001.dat", root_dir);

    inputfile = fopen(filepath, "rb");
    fread(densityArray, sizeof(densityArray[0]), sizeof(densityArray)/sizeof(densityArray[0]), inputfile);
    fclose(inputfile);

    printf("\n Total number of processes: %d", process_number);
    printf("\n local n0: %d\n", local_n0);

    for(j=0; j<alloc_local; j++) temp[j] = densityArray[j];
  }

    for(k=0; k<local_n0*n1*n2; k=k+1000){    
      if(process_rank == 0){    
        for(j=1; j<process_number; j++)  MPI_Send(&densityArray[j*local_n0*n1*n2 + k], 1000, MPI_FLOAT, j, j, MPI_COMM_WORLD);
      }
 	
      if(process_rank != 0) MPI_Recv(&temp[k], 1000, MPI_FLOAT, 0, process_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for(j=0; j<alloc_local;j++) in[j][0] = (double) temp[j];

    p = fftw_mpi_plan_dft_3d(n0, n1, n2, in, in, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    for(j=0; j<alloc_local; j++)  temp[j]         = (float) in[j][0]*pow(n0*n1*n2, -1.0);

    for(j=0; j<50; j++) printf("\n %f", temp[j]);

    for(k=0; k<local_n0*n1*n2; k=k+1000){
      if(process_rank != 0)  MPI_Send(&temp[k], 1000, MPI_FLOAT, 0, process_rank, MPI_COMM_WORLD);

      if(process_rank ==0){        
	for(j=1; j<process_number; j++)  MPI_Recv(&Array[j*local_n0*8*lenPositionArray*lenPositionArray + k + 4], 1000, MPI_FLOAT, j, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }

    for(j=0; j<alloc_local; j++)  temp[j]         = (float) in[j][1]*pow(1.0*n0*n1*n2, -1.0);

    for(k=0; k<local_n0*n1*n2; k=k+1000){
      if(process_rank != 0)  MPI_Send(&temp[k], 1000, MPI_FLOAT, 0, process_rank, MPI_COMM_WORLD);
      
      if(process_rank ==0){
        for(j=1; j<process_number; j++)  MPI_Recv(&Array[j*local_n0*8*lenPositionArray*lenPositionArray + k + 5], 1000, MPI_FLOAT, j, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } 
    }

    fftw_free(in);  

   if(process_rank == 0){
     for(s=0; s<n0; s++){
        for(q=0; q<n1; q++){
           for(r=0; r<n2; r++){
                    k_x = kInterval*r;
                    k_y = kInterval*q;
                    k_z = kInterval*s;

                    if(k_x > NyquistWavenumber)  k_x -= lenPositionArray*kInterval;
                    if(k_y > NyquistWavenumber)  k_y -= lenPositionArray*kInterval;
                    if(k_z > NyquistWavenumber)  k_z -= lenPositionArray*kInterval;
            
                    ArrayIndex        = (int) s*n1*n2 + q*n2 + r;
 
                    kSq        = (float) pow(k_x,2) + pow(k_y,2) + pow(k_z, 2);

                    Array[s*8*n1*n2+q*8*n2+r*8 + 7]           = exp(-1.*kSq*0.5*(pow(1.,2)));

                    Array[s*8*lenPositionArray*lenPositionArray+q*8*lenPositionArray+r*8+6]           = 1.;

                    if(k_x != 0.){
	                    Array[s*8*lenPositionArray*lenPositionArray+q*8*lenPositionArray+r*8+6]       *= sin(pi*k_x*kInterval*0.5/NyquistWavenumber)/(pi*k_x*kInterval*0.5/NyquistWavenumber);}
                    if(k_y != 0.){
                        Array[s*8*lenPositionArray*lenPositionArray+q*8*lenPositionArray+r*8+6]       *= sin(pi*k_y*kInterval*0.5/NyquistWavenumber)/(pi*k_y*kInterval*0.5/NyquistWavenumber);}
                    if(k_z != 0.){
                        Array[s*8*lenPositionArray*lenPositionArray+q*8*lenPositionArray+r*8+6]       *= sin(pi*k_z*kInterval*0.5/NyquistWavenumber)/(pi*k_z*kInterval*0.5/NyquistWavenumber);}

                    Array[s*8*lenPositionArray*lenPositionArray+q*8*lenPositionArray+r*8+3]           = pow(kSq, 0.5);

                    Array[s*8*lenPositionArray*lenPositionArray+q*8*lenPositionArray+r*8+2]           = k_z;
                    Array[s*8*lenPositionArray*lenPositionArray+q*8*lenPositionArray+r*8+1]           = k_y;
                    Array[s*8*lenPositionArray*lenPositionArray+q*8*lenPositionArray+r*8+0]           = k_x;
                }
            }
        }

      output = fopen("/disk1/mjw/ZADE_MockRun/FFT_NGP_kvec_Modk_Hk_Weight_3D_zadeMocks.dat", "wb");

      struct floatOutputArray k_xArray, k_yArray, k_zArray, mod_kArray, WindowFuncArray, GaussianFilterArray, H_kReal, H_kImaginary;

      for(j=0;j<n0*n1*n2;j++){
        k_xArray.floatArray[j]            = Array[8*j+0];
        k_yArray.floatArray[j]            = Array[8*j+1];
        k_zArray.floatArray[j]            = Array[8*j+2];
        mod_kArray.floatArray[j]          = Array[8*j+3];

        H_kReal.floatArray[j]             = Array[8*j+4];
        H_kImaginary.floatArray[j]        = Array[8*j+5];

        WindowFuncArray.floatArray[j]     = Array[8*j+6]; 
        GaussianFilterArray.floatArray[j] = Array[8*j+7];
      }

      fclose(output);  

    fwrite(&k_xArray,                   sizeof(struct floatOutputArray), 1, output);
    fwrite(&k_yArray,                   sizeof(struct floatOutputArray), 1, output);
    fwrite(&k_zArray,                   sizeof(struct floatOutputArray), 1, output);
    fwrite(&mod_kArray,                 sizeof(struct floatOutputArray), 1, output);
    fwrite(&H_kReal,                    sizeof(struct floatOutputArray), 1, output);
    fwrite(&H_kImaginary,               sizeof(struct floatOutputArray), 1, output);
    fwrite(&WindowFuncArray,            sizeof(struct floatOutputArray), 1, output);
    fwrite(&GaussianFilterArray,        sizeof(struct floatOutputArray), 1, output);
  }

    MPI_Finalize();
}
