#PBS -S /bin/bash                                                                                              
#PBS -N chi2_run
#PBS -V
#PBS -p 1023
#PBS -l nodes=1:ppn=1                                                                                     
#PBS -l walltime=00:30:00
#PBS -l mem=5MB

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
  export outputdir=/home/mjw/HOD_MockRun/W1_Spectro_V7_7.1
  export mask_Qldir=/home/mjw/HOD_MockRun/W1_Spectro_V7_8
  export LOZ=0.7
  export HIZ=1.2
  export FIELDFLAG=1
  export d0=1000

  rm -r /home/mjw/IO_lock/

  gcc -std=gnu11 -o chi2.o Scripts/driver_likelihood.c -lfftw3 -lm -lgsl -lgslcblas -fopenmp
}

test

DIR="$HOME/HOD_MockRun/Scripts/"
cd $DIR

export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch
export GSL_RNG_SEED=123
export GSL_RNG_TYPE=taus

export OMP_NUM_THREADS=1 # Threads = allocated processors.
cd .. 

for k in 0.6
#for k in $(seq 0.2 0.2 0.8)
  do
    echo "Analysing k_max of $k"  

    set_lock
    
    export FILE=$outputdir"/chi2_log/chi2_d0_"$d0"_W"$FIELDFLAG"_"$LOZ"_"$HIZ"_kmax_"$k".log"
      
    ./chi2.o $d0 $FIELDFLAG $LOZ $HIZ $k #> $FILE 2>&1

    ## rm -r /home/mjw/IO_lock/ ## For testing. 
    
    #if [$? -neq 0]
    #then
    #    rm -r /home/mjw/IO_lock/
    #fi
    
    #if [[ $(tr -d "\r\n" < $FILE | wc -c) -eq 0 ]]; then 
    #  printf "\n%s" "$FILE" >> $outputdir/chi2_log/chi2_stderr.pbs
    #else
    #    date >> $FILE
    #fi
done
