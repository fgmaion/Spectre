tot_np=0
fre_np=0
use_np=0

string_proccount(){
  IFS=", " read -ra numArray <<<$1

  unset count

  for word in "${numArray[@]}";
    do
      (( ${#word} == 1 )) && ((++count)) || count=$(( count +  ${word#*-} - ${word%-*} + 1))
    done

  use_np=`expr $use_np + $count`
}

for i in 30 31 32 46 47 48
  do
    jobs="$(pbsnodes worker0$i | grep 'jobs' | grep -v 'status')" 
    xx="${jobs##*= }"

    np="$(pbsnodes worker0$i | grep 'np')"
    np="${np##*= }"

    if [[ -n "${xx// }" ]]; then # check that jobs is not empty. 
        xx=$(echo $xx | sed -e 's:\/.\.head\.cluster\.loc::g' | sed -e 's:\/..\.head\.cluster\.loc::g' | sed -e 's:\/...\.head\.cluster\.loc::g' | sed -e 's:\/....\.head\.cluster\.loc::g' | sed -e 's:\/.\[.\]\.head\.cluster\.loc::g' | sed -e 's:\/..\[.\]\.head\.cluster\.loc::g' | sed -e 's:\/...\[.\]\.head\.cluster\.loc::g' | sed -e 's:\/....\[.\]\.head\.cluster\.loc::g' | sed -e 's:\/.\[..\]\.head\.cluster\.loc::g' | sed -e 's:\/..\[..\]\.head\.cluster\.loc::g' | sed -e 's:\/...\[..\]\.head\.cluster\.loc::g' | sed -e 's:\/....\[..\]\.head\.cluster\.loc::g')
       
      string_proccount $xx
    fi
      
    echo $np, $xx
    
    tot_np=`expr $tot_np + $np`
    fre_np=`expr $tot_np - $use_np`
done

echo
echo "Total processors:" $tot_np, "Free processors:" $fre_np "(Used: "$use_np")"

njobs="$(echo "scale=2; $fre_np / 8" | bc)"
njobs=$(echo $njobs | perl -nl -MPOSIX -e 'print floor($_);')

echo "Optimum number of jobs: "$njobs "(for 8 proc per job)"

nmocks=153
nruns=4 # two fields, two z-slices. 

njobs_perrun=$(python -c "import math; print int(math.floor($njobs/4))")

echo "Jobs per run: "$njobs_perrun

export nmocks_perjob=$(python -c "import math; print int(math.ceil($nmocks/$njobs_perrun))") 
echo "Mocks per job: "$nmocks_perjob

rm -rf /home/mjw/IO_lock

rm -f /home/mjw/HOD_MockRun/chi2_log/*.log
rm -f /home/mjw/HOD_MockRun/chi2_log/*.pbs
rm -f /home/mjw/HOD_MockRun/chi2_log/*.sub

rm -f /home/mjw/HOD_MockRun/pk_log/*.log
rm -f /home/mjw/HOD_MockRun/pk_log/*.pbs
rm -f /home/mjw/HOD_MockRun/pk_log/*.sub

rm -f /home/mjw/HOD_MockRun/W1_Spectro_V7_4/mocks_v1.7/dc_shifts/W4/mock_0.8_1.0_512.dat

touch /home/mjw/HOD_MockRun/chi2_log/chi2.sub
touch /home/mjw/HOD_MockRun/pk_log/pk.sub

for k in 0.6 0.8
  do         
    export LOZ=$k
    export HIZ=$(python -c "print $LOZ+0.2")
    
    for j in 1 4
      do
        export FIELDFLAG=$j

        for i in 1000 10 6 4
          do
            export d0=$i  
            
            job_id="$(qsub chi2.sh)"
            
            printf "%.1lf \t %.1lf \t W%d \t d0=% 5d, chi2_run:  %s \n" "$LOZ" "$HIZ" "$FIELDFLAG" "$d0" "$job_id" | tee -a /home/mjw/HOD_MockRun/chi2_log/chi2.sub
        done
            
        #for i in $(seq 0 $((njobs_perrun)))
        #  do
        #      export mock_start=$((1 + i*nmocks_perjob))
        #      export mock_end=$((mock_start + njobs_perrun))
        #
        #      job_id="$(qsub pk.sh)"
        #
        #      printf "%.1lf \t %.1lf \t W%d \t %03d \t pk_run:  %s \n" "$LOZ" "$HIZ" "$FIELDFLAG" "$mock_start" "$job_id" | tee -a /home/mjw/HOD_MockRun/pk_log/pk.sub
        #done
    done
done
