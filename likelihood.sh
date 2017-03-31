#PBS -S /bin/bash                                                                                              
#PBS -N likelihood.sh
#PBS -l nodes=1:ppn=1                                                                                      
#PBS -l walltime=06:00:00
#PBS -e /home/mjw/HOD_MockRun/chi2_stderr.pbs                                                                       
#PBS -o /home/mjw/HOD_MockRun/chi2_stdout.pbs                                                                                                     

DIR="$HOME/HOD_MockRun/Scripts/"
cd $DIR

export d0=10
export LOZ=0.8
export HIZ=1.0
export FIELDFLAG=1

#export LOZ=0.8
#export HIZ=1.0
#export FIELDFLAG=1

#export LOZ=0.6
#export HIZ=0.8
#export FIELDFLAG=4

#export LOZ=0.8
#export HIZ=1.0
#export FIELDFLAG=4

export BRANCH=$(git symbolic-ref --short HEAD) # current Git branch

export GSL_RNG_SEED=123
export GSL_RNG_TYPE=taus

export OMP_NUM_THREADS=1 # Threads = processors.

cd .. 

rm /home/mjw/HOD_MockRun/likelihood_d0_$d0"_W"$FIELDFLAG"_"$LOZ"_"$HIZ.log || true # or true. 

# -g: gnu debug; -w: no warnings; -o2/-03: optimization level; -DHCUBATURE; Scripts/cubature/hcubature.c; SPRNG: -lsprng -lgmp
# -std=gnu99 (for C99 with GNU extensions; https://gcc.gnu.org/onlinedocs/gcc-5.1.0/gcc/Standards.html);
# current default standard is equivalent to -std=gnu90, which is the 1989/1990 standard with GNU-specific extensions. gcc 5.1.0 (2015-04-22) changed from gnu90 to gnu11.
# -L/home/mjw/gperftools-2.5/lib -ltcmalloc -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
gcc -std=gnu11 -w -O2 -o likelihood.o Scripts/driver_likelihood.c -lfftw3 -lm -lgsl -lgslcblas

#./likelihood.o $d0 $FIELDFLAG $LOZ $HIZ 0.1 >> /home/mjw/HOD_MockRun/likelihood_d0_$d0"_W"$FIELDFLAG"_"$LOZ"_"$HIZ.log 2>&1
./likelihood.o $d0 $FIELDFLAG $LOZ $HIZ 0.2 >> /home/mjw/HOD_MockRun/likelihood_d0_$d0"_W"$FIELDFLAG"_"$LOZ"_"$HIZ.log 2>&1
./likelihood.o $d0 $FIELDFLAG $LOZ $HIZ 0.4 >> /home/mjw/HOD_MockRun/likelihood_d0_$d0"_W"$FIELDFLAG"_"$LOZ"_"$HIZ.log 2>&1
./likelihood.o $d0 $FIELDFLAG $LOZ $HIZ 0.6 >> /home/mjw/HOD_MockRun/likelihood_d0_$d0"_W"$FIELDFLAG"_"$LOZ"_"$HIZ.log 2>&1
./likelihood.o $d0 $FIELDFLAG $LOZ $HIZ 0.8 >> /home/mjw/HOD_MockRun/likelihood_d0_$d0"_W"$FIELDFLAG"_"$LOZ"_"$HIZ.log 2>&1
