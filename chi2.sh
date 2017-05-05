#PBS -S /bin/bash                                                                                              
#PBS -N chi2_run
#PBS -V
#PBS -p 1023
#PBS -l nodes=1:ppn=2                                                                                     
#PBS -l walltime=00:90:00
#PBS -l mem=5MB
#PBS -e /home/mjw/HOD_MockRun/W1_Spectro_V7_7/chi2_log/chi2_stderr.pbs                                                                       
#PBS -o /home/mjw/HOD_MockRun/W1_Spectro_V7_7/chi2_log/chi2_stdout.pbs                                                                                                     

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
  export outputdir=/home/mjw/HOD_MockRun/W1_Spectro_V7_7
  export LOZ=0.6
  export HIZ=0.9
  export FIELDFLAG=1
  export d0=1000

  rm -r /home/mjw/IO_lock/
}

test

DIR="$HOME/HOD_MockRun/Scripts/"
cd $DIR

export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch
export GSL_RNG_SEED=123
export GSL_RNG_TYPE=taus

export OMP_NUM_THREADS=2 # Threads = allocated processors.
cd .. 

# -g: gnu debug; -w: no warnings; -o2/-03: optimization level; -DHCUBATURE; Scripts/cubature/hcubature.c; SPRNG: -lsprng -lgmp
# -std=gnu99 (for C99 with GNU extensions; https://gcc.gnu.org/onlinedocs/gcc-5.1.0/gcc/Standards.html);
# current default standard is equivalent to -std=gnu90, which is the 1989/1990 standard with GNU-specific extensions. gcc 5.1.0 (2015-04-22) changed from gnu90 to gnu11.
# -L/home/mjw/gperftools-2.5/lib -ltcmalloc -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
gcc -std=gnu11 -w -O2 -o chi2.o Scripts/driver_likelihood.c -lfftw3 -lm -lgsl -lgslcblas -fopenmp

for k in $(seq 0.2 0.2 0.8)
  do
    echo "Analysing k_max of $k"  

    set_lock
    
    export FILE=$outputdir"/chi2_log/chi2_d0_"$d0"_W"$FIELDFLAG"_"$LOZ"_"$HIZ"_kmax_"$k".log"
      
    ./chi2.o $d0 $FIELDFLAG $LOZ $HIZ $k # > $FILE 2>&1

    if [[ $(tr -d "\r\n" < $FILE | wc -c) -eq 0 ]]; then 
      printf "\n%s" "$FILE" >> $outputdir/chi2_log/chi2_stderr.pbs
    else
        date >> $FILE
    fi
done
