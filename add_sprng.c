#include <mpi.h>

#include "/home/mjw/sprng2.0/gsl-sprng.h" // gsl-sprng interface: https://darrenjw.wordpress.com/tag/gsl/

int test(){
  int i, k, po;

  gsl_rng *r;

  // MPI_Comm_rank(MPI_COMM_WORLD, &k);

  r = gsl_rng_alloc(gsl_rng_sprng20);

  for (i=0; i<10; i++){
    po = gsl_ran_poisson(r, 2.0);

    printf("Process %d, random number %d: %d\n", k, i+1, po);
  }

  return 0;
}
