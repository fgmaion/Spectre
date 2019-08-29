#PBS -S /bin/bash
#PBS -N pk_run
#PBS -V 
#PBS -l nodes=1:ppn=1
#PBS -l walltime=0:120:00
#PBS -l mem=5MB
#PBS -p 1023   

set_lock(){
    locked=1

    while [ $locked = 1 ]
    do
        sleep  1

        locked=0

        DIR=/home/mjw/IO_lock

        if mkdir $DIR; then
            echo "Locking succeeded"  > $DIR/lock.dat
        else
            locked=1
        fi
    done
}

DIR="/global/homes/m/mjwilson/Spectrum-MC"
cd $DIR

export OMP_NUM_THREADS=1                       # Threads = processors.
export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch

export GSL_RNG_TYPE="taus"
export GSL_RNG_SEED=123

# -pedantic -std=gnu11 -Wall -Wextra
gcc -o pk.o driver_pk.c -I/opt/cray/pe/fftw/3.3.8.2/haswell/include -L/opt/cray/pe/fftw/3.3.8.2/haswell/lib -fopenmp -lfftw3_omp -lfftw3 -lm -I/global/common/sw/cray/cnl7/haswell/gsl/2.5/intel/19.0.3.199/7twqxxq/include -L/global/common/sw/cray/cnl7/haswell/gsl/2.5/intel/19.0.3.199/7twqxxq/lib -lgsl -lgslcblas 

./pk.o 1 0.6 0.9
