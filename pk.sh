#PBS -S /bin/bash
#PBS -N pk_run
#PBS -V 
#PBS -l nodes=1:ppn=2
#PBS -l walltime=0:25:00
#PBS -l mem=5MB
#PBS -p 1023   
#PBS -e /home/mjw/HOD_MockRun/W1_Spectro_V7_7/pk_log/pk_stderr.pbs
#PBS -o /home/mjw/HOD_MockRun/W1_Spectro_V7_7/pk_log/pk_stdout.pbs

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

DIR="/home/mjw/HOD_MockRun/Scripts/"
cd $DIR

export OMP_NUM_THREADS=2 # Threads = processors.
export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch

export GSL_RNG_TYPE="taus"
export GSL_RNG_SEED=123


test(){
  export outputdir=/home/mjw/HOD_MockRun/W1_Spectro_V7_7

  export mock_start=153   # not zero!
  export nmocks_perjob=1 
  export LOZ=0.9
  export HIZ=1.2
  export FIELDFLAG=1
}

## test
  
cd /home/mjw/HOD_MockRun/

set_lock

date >> $outputdir/pk_log/"pk_W"$FIELDFLAG"_"$LOZ"_"$HIZ"_"$mock_start".log"

# -g: gnu debug; -w: no warnings; -o2/-03: optimization level; -DHCUBATURE; Scripts/cubature/hcubature.c;
# -std=gnu99 (for C99 with GNU extensions; https://gcc.gnu.org/onlinedocs/gcc-5.1.0/gcc/Standards.html);
# current default standard is equivalent to -std=gnu90, which is the 1989/1990 standard with GNU-specific
# extensions. gcc 5.1.0 (2015-04-22) changed from gnu90 to gnu11.

# -O2: speed, -g: debug; -Werror
gcc -O2 -Wall -pedantic -Wextra -std=gnu11 -o pk.o Scripts/driver_pk_d0.c -fopenmp -lfftw3_omp -lfftw3 -lm  -lgsl -lgslcblas

# /home/ert/local/bin/valgrind --tool=memcheck --leak-check=full /home/mjw/HOD_MockRun/pk.o $FIELDFLAG $LOZ $HIZ $mock_start $nmocks_perjob
# /home/ert/local/bin/valgrind --tool=massif --stacks=yes

./pk.o $FIELDFLAG $LOZ $HIZ $mock_start $nmocks_perjob >> $outputdir/pk_log/"pk_W"$FIELDFLAG"_"$LOZ"_"$HIZ"_"$mock_start".log" 2>&1
