#PBS -S /bin/bash
#PBS -N xiq_run
#PBS -V
#PBS -p 1023
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00                                                                                      

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
    export version=9
    export oldversion="/home/mjw/HOD_MockRun/W1_Spectro_V7_"`expr $version - 1`

    export outputdir="/home/mjw/HOD_MockRun/W1_Spectro_V7_"$version
    export mask_Qldir=/home/mjw/HOD_MockRun/W1_Spectro_V7_2
     
    ## Interactive run with: qsub -I -o $outputdir/chi2_log/chi2_stdout.pbs -e $outputdir/chi2_log/chi2_stderr.pbs chi2.sh
    export LOZ=0.6
    export HIZ=0.8
    export FIELDFLAG=1
    export d0=1000
    export RES=2
    export NPROCESSOR=4
    
    MAX_LOGS=(2. 20. 2000. 2. 20. 4000.)
    MAX_LOG=${MAX_LOGS[$RES]}

    ## NLOGBINS IN R x NLINBINS IN MU (20).
    export NBINS=$(python -c "from math import ceil, log10; print int(40*ceil((log10(4000.) - log10(0.001))/log10(1.05)))")

    ## -I /share/star/include/ -L/share/star/lib/
    gcc -o xiq_$RES.o Scripts/driver_Qmultipoles.c -fopenmp -lfftw3 -lm  -lgsl -lgslcblas -lcfitsio -D "NBINS=${NBINS}"
    
    rm -rf /home/mjw/IO_lock/
}

## test 

DIR="$HOME/HOD_MockRun/Scripts"
cd $DIR

export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch
export GSL_RNG_SEED=123
export GSL_RNG_TYPE=taus

export OMP_NUM_THREADS=$NPROCESSOR # Threads = allocated processors.
cd ..

set_lock

export FILE=$outputdir"/xiq_log/xiq_W"$FIELDFLAG"_"$LOZ"_"$HIZ"_"$RES".log"

./xiq_$RES.o $FIELDFLAG $LOZ $HIZ $RES > $FILE 2>&1

if [[ $(tr -d "\r\n" < $FILE | wc -c) -eq 0 ]]; then
    printf "\n%s" "$FILE" >> $outputdir/xiq_log/xiq_stderr.pbs
else
    date >> $FILE
fi
