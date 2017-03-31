#PBS -S /bin/bash
#PBS -V ## All nodes get correct environment.
#PBS -l nodes=1:ppn=16
#PBS -l walltime=0:100:00

#PBS -e /home/mjw/HOD_MockRun/pk_stderr.pbs
#PBS -o /home/mjw/HOD_MockRun/pk_stdout.pbs

DIR="/home/mjw/HOD_MockRun/Scripts/"
cd $DIR

export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch

export GSL_RNG_TYPE="taus"
export GSL_RNG_SEED=123

export LOZ=0.6
export HIZ=0.8
export FIELDFLAG=4

cd ..

export OMP_NUM_THREADS=16 # Threads = processors.

rm -f "pk_W"$FIELDFLAG"_"$LOZ"_"$HIZ".log"

date >> "pk_W"$FIELDFLAG"_"$LOZ"_"$HIZ".log"

# -g: gnu debug; -w: no warnings; -o2/-03: optimization level; -DHCUBATURE; Scripts/cubature/hcubature.c;
# -std=gnu99 (for C99 with GNU extensions; https://gcc.gnu.org/onlinedocs/gcc-5.1.0/gcc/Standards.html);
# current default standard is equivalent to -std=gnu90, which is the 1989/1990 standard with GNU-specific extensions. gcc 5.1.0 (2015-04-22) changed from gnu90 to gnu11.
gcc -std=gnu11 -w -O2 -o pk.o Scripts/driver_pk_d0.c -fopenmp -lfftw3_omp -lfftw3 -lm  -lgsl -lgslcblas

./pk.o $FIELDFLAG $LOZ $HIZ  >> "pk_W"$FIELDFLAG"_"$LOZ"_"$HIZ".log" 2>&1
