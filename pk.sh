#PBS -S /bin/bash                                                                                              
#PBS -l nodes=3:ppn=16
#PBS -l walltime=0:10:00
#PBS -e /home/mjw/HOD_MockRun/pk_errors.pbs                                                                       
#PBS -o /home/mjw/HOD_MockRun/pk_messages.pbs                                                                                                     

DIR="/home/mjw/HOD_MockRun/Scripts/"
cd $DIR

branch=$(git symbolic-ref --short HEAD) # current Git branch

echo "GIT BRANCH:" $branch

cd .. 

export OMP_NUM_THREADS=16 # Threads = processors.

# -g: gnu debug; -w: no warnings; -o2/-03: optimization level; -DHCUBATURE; Scripts/cubature/hcubature.c; SPRNG: -lsprng -lgmp 
# -std=gnu99 (for C99 with GNU extensions; https://gcc.gnu.org/onlinedocs/gcc-5.1.0/gcc/Standards.html);
# current default standard is equivalent to -std=gnu90, which is the 1989/1990 standard with GNU-specific extensions. gcc 5.1.0 (2015-04-22) changed from gnu90 to gnu11. 
gcc -std=gnu11 -g -w -O2 -o driver_pk_d0.o Scripts/driver_pk_d0.c -fopenmp -lfftw3_omp -lfftw3 -lm -lgsl -lgslcblas -L/home/mjw/gperftools-2.5/lib -ltcmalloc -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free

GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_pk_d0.o 1 0.6 0.9
#gdb driver_pk_d0.o
