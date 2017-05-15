#PBS -S /bin/bash
#PBS -N xiq_run
#PBS -V
#PBS -p 1023
#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00                                                                                      

set_lock(){
    locked=1

    while [ $locked = 1 ]
    do
        sleep  1

        locked=0

        DIR=/home/mjw/IO_lock

        if mkdir $DIR; then
            echo "New lock placed"
        else
            echo "Still locked"

            locked=1
        fi
    done
}

test(){
    ## Interactive run with: qsub -I -o $outputdir/chi2_log/chi2_stdout.pbs -e $outputdir/chi2_log/chi2_stderr.pbs chi2.sh
    export outputdir=/home/mjw/HOD_MockRun/W1_Spectro_V7_7
    export mask_Qldir=/home/mjw/HOD_MockRun/W1_Spectro_V7_2
    export LOZ=0.6
    export HIZ=0.9
    export FIELDFLAG=1
    export d0=1000
    export NBINS=15280
    export RES=0
    
    gcc -w -O2 -fopenmp -o xi_q.o Scripts/driver_Qmultipoles.c -lfftw3 -lm  -lgsl -lgslcblas -D "NBINS=${NBINS}"
    
    rm -r /home/mjw/IO_lock/
}

## test 

DIR="$HOME/HOD_MockRun/Scripts"
cd $DIR

export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch
export GSL_RNG_SEED=123
export GSL_RNG_TYPE=taus

export OMP_NUM_THREADS=8 # Threads = allocated processors.
cd ..

set_lock

export FILE=$outputdir"/xiq_log/xiq_W"$FIELDFLAG"_"$LOZ"_"$HIZ"_"$RES".log"

./xi_q.o $FIELDFLAG $LOZ $HIZ $RES # > $FILE 2>&1

if [[ $(tr -d "\r\n" < $FILE | wc -c) -eq 0 ]]; then
    printf "\n%s" "$FILE" >> $outputdir/xiq_log/xiq_stderr.pbs
else
    date >> $FILE
fi
