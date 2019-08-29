#PBS -S /bin/bash
#PBS -N Fisher_run
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -l walltime=0:20:00
#PBS -p 1023
#PBS -e /home/mjw/HOD_MockRun/W1_Spectro_V7_7/fisher/Fisher_stderr.pbs
#PBS -o /home/mjw/HOD_MockRun/W1_Spectro_V7_7/fisher/Fisher_stdout.pbs 

test(){
    export outputdir=/home/mjw/HOD_MockRun/W1_Spectro_V7_7.1
    export mask_Qldir=/home/mjw/HOD_MockRun/W1_Spectro_V7_2
    export LOZ=0.6
    export HIZ=0.9
    export FIELDFLAG=4
    export d0=1000
    export KMAX=0.8
    
    rm -r /home/mjw/IO_lock/
}

test 

DIR="$HOME/HOD_MockRun/Scripts/"
cd $DIR

export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch
export GSL_RNG_SEED=123
export GSL_RNG_TYPE=taus

export OMP_NUM_THREADS=1 # Threads = allocated processors.
cd .. 

# -g: gnu debug; -w: no warnings; -o2/-03: optimization level; -DHCUBATURE; Scripts/cubature/hcubature.c;
# -std=gnu99 (for C99 with GNU extensions; https://gcc.gnu.org/onlinedocs/gcc-5.1.0/gcc/Standards.html);
# current default standard is equivalent to -std=gnu90, which is the 1989/1990 standard with GNU-specific
# extensions. gcc 5.1.0 (2015-04-22) changed from gnu90 to gnu11.
gcc -pedantic -std=gnu11 -O2 -o Fisher.o Scripts/driver_Fisher.c -fopenmp -lfftw3_omp -lfftw3 -lm  -lgsl -lgslcblas

./Fisher.o $d0 $FIELDFLAG $LOZ $HIZ $KMAX


