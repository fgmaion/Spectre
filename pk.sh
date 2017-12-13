#PBS -S /bin/bash
#PBS -N pk_run
#PBS -V 
#PBS -l nodes=1:ppn=4
#PBS -l walltime=0:100:00
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

DIR="/home/mjw/HOD_MockRun/Scripts/"
cd $DIR

export OMP_NUM_THREADS=4 # Threads = processors.
export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch

export GSL_RNG_TYPE="taus"
export GSL_RNG_SEED=123

cd ..

test(){
  ## Interactive run with: qsub -I -o $outputdir/pk_log/pk_stdout.pbs -e $outputdir/pk_log/pk_stderr.pbs pk.sh
  export outputdir=/home/mjw/HOD_MockRun/W1_Spectro_V7_9/
  export mock_start=1   # not zero!
  export nmocks_perjob=1 
  export LOZ=0.6
  export HIZ=0.9
  export FIELDFLAG=1

  rm -rf /home/mjw/IO_lock/

  gcc -Wall -pedantic -Wextra -std=gnu11 -o pk.o Scripts/driver_pk_d0.c -fopenmp -lfftw3_omp -lfftw3 -lm  -lgsl -lgslcblas
}

test

set_lock

date >> $outputdir/pk_log/"pk_W"$FIELDFLAG"_"$LOZ"_"$HIZ"_"$mock_start".log"

# /home/ert/local/bin/valgrind --tool=memcheck --leak-check=full /home/mjw/HOD_MockRun/pk.o $FIELDFLAG $LOZ $HIZ $mock_start $nmocks_perjob
# /home/ert/local/bin/valgrind --tool=massif --stacks=yes

./pk.o $FIELDFLAG $LOZ $HIZ $mock_start $nmocks_perjob  # >> $outputdir/pk_log/"pk_W"$FIELDFLAG"_"$LOZ"_"$HIZ"_"$mock_start".log" 2>&1

#if [$? -neq 0]
#then
#    rm -r /home/mjw/IO_lock/
#fi
