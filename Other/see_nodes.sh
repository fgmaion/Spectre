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

