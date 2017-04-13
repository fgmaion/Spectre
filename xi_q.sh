#PBS -S /bin/bash
#PBS -N xiq_run
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00                                                                                      
#PBS -e /home/mjw/HOD_MockRun/xiq_log/xiq_stderr.pbs                                                                       
#PBS -o /home/mjw/HOD_MockRun/xiq_log/xiq_stdout.pbs
                                                                                                     
DIR="$HOME/HOD_MockRun/"
cd $DIR

# 0: hihi res. 
# sampling_frac  =                        0.7;                                                                                                                                               # maxlog         =               log10(  2.0);                                                                                                                                                                  
# 1: hi res.
# sampling_frac  =                        0.1;                                                                                                                                               # maxlog         =               log10( 20.0);                                                                                                                                                                  
# 2: lo res.
# sampling_frac  =                       0.01;                                                                                                                                               # maxlog         =              log10(2000.0);  

gcc -w -O2 -fopenmp -o xi_q.o Scripts/driver_Qmultipoles.c -lfftw3 -lm  -lgsl -lgslcblas

## 0.6 < z < 0.9;    Last argument is the crucial one, see paircount.c for an explanation. 
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.9 1.0 2.0 1.0 0      # W1, 0.6-0.9, hihiRes
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.9 1.0 20.0 0.1 1     # W1, 0.6-0.9,   hiRes 
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.9 1.0 2000.0 0.01 2  # W1, 0.6-0.9,   loRes 

## 0.9 < z < 1.2
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.9 1.2 1.0 2.0  1.0 0      # W1, 0.9-1.2, hihiRes
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.9 1.2 1.0 20.0 0.1 1      # W1, 0.9-1.2,   hiRes
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.9 1.2 1.0 2000.0 0.01 2   # W1, 0.9-1.2,   loRes


## 0.6 < z < 0.8
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.8 1.0 2.0  1.0 0      # W1, 0.8-1.0, hihiRes                                                                     
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.8 1.0 20.0 0.1 1      # W1, 0.8-1.0,   hiRes                                                                     
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.8 1.0 2000.0 0.01 2   # W1, 0.8-1.0,   loRes

## 0.8 < z < 1.0
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.8 1.0 1.0 2.0  1.0 0      # W1, 0.8-1.0, hihiRes                                       
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.8 1.0 1.0 20.0 0.1 1      # W1, 0.8-1.0,   hiRes                                      
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.8 1.0 1.0 2000.0 0.01 2   # W1, 0.8-1.0,   loRes  


## W4 field
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 4 1000.0 0.6 0.9 1.0 2.0  1.0 0      # W4, 0.6-0.9, hihiRes
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 4 1000.0 0.6 0.9 1.0 20.0 0.1 1      # W4, 0.6-0.9,   hiRes
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 4 1000.0 0.6 0.9 1.0 2000.0 0.01 2   # W4, 0.6-0.9,   loRes

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 4 1000.0 0.9 1.2 1.0 2.0 1.0 0       # W4, 0.9-1.2, hihiRes
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 4 1000.0 0.9 1.2 1.0 20.0 0.1 1      # W4, 0.9-1.2,   hiRes
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 4 1000.0 0.9 1.2 1.0 2000.0 0.01 2   # W4, 0.9-1.2,   loRes

GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./xi_q.o 4 0.6 0.8    2.0 1.00 0   # W4, 0.6-0.8, hihiRes
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./xi_q.o 4 0.6 0.8   20.0 0.10 1   # W4, 0.6-0.8,   hiRes
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./xi_q.o 4 0.6 0.8 2000.0 0.01 2   # W4, 0.6-0.8,   loRes

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./xi_q.o 4 0.8 1.0    2.0 1.00 0   # W4, 0.8-1.0, hihiRes
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./xi_q.o 4 0.8 1.0   20.0 0.10 1   # W4, 0.8-1.0,   hiRes
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./xi_q.o 4 0.8 1.0 2000.0 0.01 2   # W4, 0.8-1.0,   loRes 

### JOINT-FIELDS. 
## Load mask hard coded in this case.  
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.9 1.0 2.0    1.000 3  # W1 W4, 0.6-0.9, hihiRes                                                               
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.9 1.0 20.0   0.004 4  # W1 W4, 0.6-0.9,   hiRes                                              
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.9 1.0 2000.0 0.005 5  # W1 W4, 0.6-0.9,   loRes                                                                  

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.9 1.2 1.0 2.0    1.000 3  # W1 W4, 0.9-1.2, hihiRes                                                               
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.9 1.2 1.0 20.0   0.004 4  # W1 W4, 0.9-1.2,   hiRes                                                               
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.9 1.2 1.0 2000.0 0.005 5  # W1 W4, 0.9-1.2,   loRes 

## 0.6 < z < 0.8
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.8 1.0 2.0    1.000 3  # W1 W4, 0.6-0.9, hihiRes                                                                  
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.8 1.0 20.0   0.004 4  # W1 W4, 0.6-0.9,   hiRes                                                                  
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.6 0.8 1.0 2000.0 0.005 5  # W1 W4, 0.6-0.9,   loRes                                                                    
## 0.8 < z < 1.0
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.8 1.0 1.0 2.0    1.000 3  # W1 W4, 0.9-1.2, hihiRes                                                                  
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.8 1.0 1.0 20.0   0.004 4  # W1 W4, 0.9-1.2,   hiRes                                                                  
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_Qmultipoles.o 1 1000.0 0.8 1.0 1.0 2000.0 0.005 5  # W1 W4, 0.9-1.2,   loRes 
