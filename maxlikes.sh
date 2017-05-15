#PBS -S /bin/bash
#PBS -N maxlikes_run
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -l walltime=0:20:00
#PBS -p 1023

test(){
  ## Interactive run with: qsub -I -o $outputdir/maxlikes/maxlikes_stdout.pbs -e $outputdir/maxlikes/maxlikes_stderr.pbs maxlikes.sh
  export outputdir=/home/mjw/HOD_MockRun/W1_Spectro_V7_7
  export mask_Qldir=/home/mjw/HOD_MockRun/W1_Spectro_V7_2
  export LOZ=0.6
  export HIZ=0.9
  export FIELDFLAG=1
  export d0=1000
  export KMAX=0.8
  
  rm -r /home/mjw/IO_lock/

  ## likelihood as shared. 
  gcc -shared -fPIC -std=gnu11 -w -O2 -o likelihood.so /home/mjw/HOD_MockRun/Scripts/driver_likelihood_ctypes.c -fopenmp -lfftw3_omp -lfftw3 -lm  -lgsl -lgslcblas
}

## test

DIR="$HOME/HOD_MockRun/Scripts/"
cd $DIR

export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch
export GSL_RNG_SEED=123
export GSL_RNG_TYPE=taus

export OMP_NUM_THREADS=1 # Threads = allocated processors.
cd .. 

python Scripts/get_maxlikes.py


