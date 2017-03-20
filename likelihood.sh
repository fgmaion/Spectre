#PBS -S /bin/bash                                                                                              
#PBS -l nodes=1:ppn=1:walltime=48:10:00                                                                                      
#PBS -e /home/mjw/HOD_MockRun/likelihood_errors.pbs                                                                       
#PBS -o /home/mjw/HOD_MockRun/likelihood_messages.pbs                                                                                                     

DIR="$HOME/HOD_MockRun/"
cd $DIR

# no warnings
gcc -w -DHCUBATURE -fopenmp -lfftw3_omp -o driver_likelihood.o Scripts/cubature/hcubature.c Scripts/driver_likelihood.c -lfftw3 -lm -I/home/ert/local/gsl/current/gcc-4.7.2/include/ -L/home/ert/local/gsl/current/gcc-4.7.2/lib/ -lgsl -lgslcblas

export OMP_NUM_THREADS=1  # Threads = processors.  

GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.6 0.9 1.0 8

## for mean multipoles
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.6 0.9 1.0 8 # > likelihood_output.txt
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.9 1.2 1.0 8 > likelihood_output.txt                    
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.6 0.9 1.0 8 > likelihood_output.txt                    
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.9 1.2 1.0 8 > likelihood_output.txt

## for k_max
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.6 0.9 1.0 1 > likelihood_output.txt                                                   #GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.9 1.2 1.0 1 > likelihood_output.txt                                                   #GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.6 0.9 1.0 1 > likelihood_output.txt                                                   #GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.9 1.2 1.0 1 > likelihood_output.txt

#for((i=1; i<26; i++))
#  do  
#    qsub  singleRun_likelihood_kmax0.4.sh -v NUM=$i

#    sleep 20 

#    qsub singleRun_likelihood_kmax0.6.sh -v NUM=$i

#    sleep 20

#    qsub singleRun_likelihood_kmax0.8.sh -v NUM=$i

#    sleep 60  
#  done

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.6 0.9 1.0 8                                                                           #GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.9 1.2 1.0 8                                                                           #GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.6 0.9 1.0 8                                                                           #GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.9 1.2 1.0 8

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.6 0.9 1.0 6
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.9 1.2 1.0 6
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.6 0.9 1.0 6
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.9 1.2 1.0 6

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.6 0.9 1.0 4
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.9 1.2 1.0 4
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.6 0.9 1.0 4
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.9 1.2 1.0 4

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 10.0 0.6 0.9 1.0 8
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 10.0 0.6 0.9 1.0 4
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 10.0 0.6 0.9 1.0 4
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 10.0 0.6 0.9 1.0 4

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 6.0 0.6 0.9 1.0 8
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 6.0 0.6 0.9 1.0 4
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 6.0 0.6 0.9 1.0 4
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 6.0 0.6 0.9 1.0 4

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 4.0 0.6 0.9 1.0 8

## 0.8 < z < 1.0
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 4.0 0.8 1.0 1.0 8

#for each kmax value. 
#for((i=2; i<8; i=i+2))
#  do
#   echo $i 
    
   #GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.6 0.9 1.0 $i                                                                                                     
   #GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 1 1000.0 0.9 1.2 1.0 $i
   #GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.6 0.9 1.0 $i  
#   GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.9 1.2 1.0 $i     
#  done  


## ******************** W4 *********************** ##
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.6 0.9 1.0 4                                                                                                                                                          
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.6 0.9 1.0 6                                                                            
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 1000.0 0.6 0.9 1.0 8 

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 10.0 0.6 0.9 1.0 4                                                                                                                                                            
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 10.0 0.6 0.9 1.0 6                                                                                                                                                            
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 10.0 0.6 0.9 1.0 8 

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 6.0 0.6 0.9 1.0 4                                                                                                                                                             
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 6.0 0.6 0.9 1.0 6                                                                                                                                                             
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 6.0 0.6 0.9 1.0 8

#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 4.0 0.6 0.9 1.0 4
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 4.0 0.6 0.9 1.0 6
#GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./driver_likelihood.o 4 4.0 0.6 0.9 1.0 8
