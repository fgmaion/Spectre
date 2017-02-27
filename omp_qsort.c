#include "assert.h"
#include "omp.h"
#include "smartomp.h"
#include "stdlib.h"
#include "stdio.h"


int max_threads;

/* How many workers can be working at same time */
int num_of_workers = 0;
/* How many active workers currently */
int num_of_active_workers = 0;

omp_lock_t * worker_lock;

int * generate_data(int n)
{
  int i=0;
  int * data = malloc(sizeof(int) * n);
  for (i = 0; i < n; i++)
    data[i] = random();

  return data;
}


void _qsort(int *data, int n, int was_spawned)
{
  if (n >  1) {

    //printf("Thread id %d %d\n", omp_get_thread_num(), n);

    int pivot = data[0];
    int l = 0;
    int r = n;

    while (1)
      {
	while (l < n && data[++l] < pivot);
	while (data[--r] > pivot);
	if (l >= r)
	  break;
	int t = data[l];
	data[l] = data[r];
	data[r] = t;
      }

    int t = data[l - 1];
    data[l - 1] = data[0];
    data[0] = t;

    int recursion_in_parallel = 0;
    omp_set_lock(worker_lock);

    /* Limit number of threads running */
    if (num_of_active_workers < num_of_workers) {
      recursion_in_parallel = 1;
      num_of_active_workers += 2;
    }
    omp_unset_lock(worker_lock);

    if (recursion_in_parallel == 1) {
      #pragma omp parallel sections
      {
	 #pragma omp section
	_qsort(data, l - 1, 1);

	 #pragma omp section
	_qsort(data + l, n - l, 1);
      }
    } else {
      _qsort(data, l - 1, 0);
      _qsort(data + l, n - l, 0);
    }
  }
  /* Decrease num */
  if (was_spawned == 1) {
    omp_set_lock(worker_lock);
    num_of_active_workers--;
    omp_unset_lock(worker_lock);
  }
}

void verify(int *data, int n)
{
  int i;
  for (i = 1; i < n; i++)
    assert(data[i] >= data[i-1]);
}


void parallel_qsort(int *data, int n) {
  int num_of_threads = max_threads;
  omp_set_num_threads(max_threads);
  smartomp_ticket_ptr ticket =
    SMARTOMP_start_timing_get_threads("quicksort",  n, &num_of_threads);
  // Set number of threads to be used for the parallel sections

  printf("====================================\n");
  printf("Start qsort with %d/%d/%d threads\n", num_of_threads, omp_get_num_threads(), max_threads);

  omp_set_num_threads(num_of_threads);
  num_of_workers = num_of_threads;
  num_of_active_workers = 1;
  _qsort(data, n, 1);

  SMARTOMP_end_timing(ticket);
}


int main(int argv, char *argc[])
{
  if (argv < 3) {
    printf("Usage: quicksort [num_of_threads] [num_of_iterations] Exiting.\n");
    exit(1);
  }
  int test_input_size[20] = {1e5, 2e5, 3e5, 5e5, 1e6, 2e6, 3e6, 5e6, 1e7, 2e7, 1e8, 2e8,5e8,8e8, 1e9,2e9, 0}; // Last must be 0
  sscanf(argc[1], "%d", &max_threads);
  printf("Maximum num of threads: %d, %s\n", max_threads, argc[1]);
  omp_set_nested(1);
  // Set seed to get always same numbers
  srandom(1);
  int ITERATIONS = 0;
  sscanf(argc[2], "%d", &ITERATIONS);

  printf("Going to run %d iterations\n", ITERATIONS);

  int iter = 0;

  /* Init lock */
  worker_lock = malloc(sizeof(omp_lock_t));
  omp_init_lock(worker_lock);

  double starttime = omp_get_wtime();

  for(iter=0; iter<ITERATIONS; iter++) {
    int input_iter = 0;
    while(test_input_size[input_iter] != 0) {
      int N = test_input_size[input_iter++];
      int * data = generate_data(N);

      double start_time = omp_get_wtime();
      parallel_qsort(data, N);
      double end_time = omp_get_wtime();
      verify(data, N);
      free(data);

      printf("Num of active workers: %d\n", num_of_active_workers);
      assert(num_of_active_workers == 0);


      printf("Sort of %d ints took %lf secs\n", N, (end_time-start_time));
    }
  }

  double endtime = omp_get_wtime();

  printf("Total time: %lf\n", (endtime-starttime));

  return 0;
}
