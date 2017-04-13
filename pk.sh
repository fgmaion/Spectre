#PBS -S /bin/bash
#PBS -N pk_run
#PBS -V 
#PBS -l nodes=1:ppn=8
#PBS -l walltime=0:20:00
#PBS -l mem=5MB
#PBS -p 1023   
#PBS -e /home/mjw/HOD_MockRun/pk_log/pk_stderr.pbs
#PBS -o /home/mjw/HOD_MockRun/pk_log/pk_stdout.pbs

DIR="/home/mjw/HOD_MockRun/Scripts/"
cd $DIR

export OMP_NUM_THREADS=8 # Threads = processors.
export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch

export GSL_RNG_TYPE="taus"
export GSL_RNG_SEED=123

test(){
  export mock_start=1   # not zero!
  export nmocks_perjob=1 
  export LOZ=0.6
  export HIZ=0.8
  export FIELDFLAG=4
}

test
  
cd /home/mjw/HOD_MockRun/

date >> /home/mjw/HOD_MockRun/pk_log/"pk_W"$FIELDFLAG"_"$LOZ"_"$HIZ"_"$mock_start".log"

# -g: gnu debug; -w: no warnings; -o2/-03: optimization level; -DHCUBATURE; Scripts/cubature/hcubature.c;
# -std=gnu99 (for C99 with GNU extensions; https://gcc.gnu.org/onlinedocs/gcc-5.1.0/gcc/Standards.html);

# current default standard is equivalent to -std=gnu90, which is the 1989/1990 standard with GNU-specific
# extensions. gcc 5.1.0 (2015-04-22) changed from gnu90 to gnu11.
gcc -std=gnu11 -w -O2 -o pk.o Scripts/driver_pk_d0.c -fopenmp -lfftw3_omp -lfftw3 -lm  -lgsl -lgslcblas

# ${PBS_JOBID%.head.cluster*}
# /home/ert/local/bin/valgrind --tool=massif --stacks=yes
./pk.o $FIELDFLAG $LOZ $HIZ $mock_start $nmocks_perjob # >> /home/mjw/HOD_MockRun/pk_log/"pk_W"$FIELDFLAG"_"$LOZ"_"$HIZ"_"$mock_start".log" 2>&1
