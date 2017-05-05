#PBS -S /bin/bash
#PBS -N emcee_run
#PBS -V
#PBS -l nodes=1:ppn=4
#PBS -l walltime=0:20:00
#PBS -p 1023
#PBS -e /home/mjw/HOD_MockRun/W1_Spectro_V7_7/emcee_log/emcee_stderr.pbs
#PBS -o /home/mjw/HOD_MockRun/W1_Spectro_V7_7/emcee_log/emcee_stdout.pbs 

test(){
  export outputdir=/home/mjw/HOD_MockRun/W1_Spectro_V7_7
  export LOZ=0.6
  export HIZ=0.9
  export FIELDFLAG=1
  export D0=1000
  export KMAX=0.8
  
  rm -r /home/mjw/IO_lock/
}

test

DIR="$HOME/HOD_MockRun/Scripts/"
cd $DIR

export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch
export GSL_RNG_SEED=123
export GSL_RNG_TYPE=taus

export OMP_NUM_THREADS=4 # Threads = allocated processors.
cd .. 

# -g: gnu debug; -w: no warnings; -o2/-03: optimization level; -DHCUBATURE; Scripts/cubature/hcubature.c;
# -std=gnu99 (for C99 with GNU extensions; https://gcc.gnu.org/onlinedocs/gcc-5.1.0/gcc/Standards.html);
# current default standard is equivalent to -std=gnu90, which is the 1989/1990 standard with GNU-specific
# extensions. gcc 5.1.0 (2015-04-22) changed from gnu90 to gnu11.

mpirun -np $OMP_NUM_THREADS /home/mlam/anaconda/bin/python Scripts/run_emcee.py


