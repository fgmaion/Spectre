#PBS -S /bin/bash
#PBS -N pk.sh
#PBS -l nodes=1:ppn=16
#PBS -l walltime=0:30:00 
#PBS -e pk_stderr.pbs
#PBS -o pk_stdout.pbs

export DIR="/home/mjw/HOD_MockRun/Scripts/"
cd $DIR

#export LOZ=0.6
#export HIZ=0.8
#export FIELDFLAG=1

#export LOZ=0.8
#export HIZ=1.0
#export FIELDFLAG=1

#export LOZ=0.6
#export HIZ=0.8
#export FIELDFLAG=4

export LOZ=0.8
export HIZ=1.0
export FIELDFLAG=4 


export USCORE=_

export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch

export GSL_RNG_SEED=123
export GSL_RNG_TYPE=taus

export OMP_NUM_THREADS=16 # Threads = processors. 

cd .. 

rm /home/mjw/HOD_MockRun/pk_W$FIELDFLAG$USCORE$LOZ$USCORE$HIZ.log

# -g: gnu debug; -w: no warnings; -o2/-03: optimization level; -DHCUBATURE; Scripts/cubature/hcubature.c; SPRNG: -lsprng -lgmp 
# -std=gnu99 (for C99 with GNU extensions; https://gcc.gnu.org/onlinedocs/gcc-5.1.0/gcc/Standards.html);
# current default standard is equivalent to -std=gnu90, which is the 1989/1990 standard with GNU-specific extensions. gcc 5.1.0 (2015-04-22) changed from gnu90 to gnu11. 
# -L/home/mjw/gperftools-2.5/lib -ltcmalloc -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
gcc -std=gnu11 -w -O2 -o pk_d0.o Scripts/driver_pk_d0.c -fopenmp -lfftw3_omp -lfftw3 -lm -lgsl -lgslcblas

./pk_d0.o $FIELDFLAG $LOZ $HIZ >> /home/mjw/HOD_MockRun/pk_W$FIELDFLAG$USCORE$LOZ$USCORE$HIZ.log 2>&1

#gdb driver_pk_d0.o
