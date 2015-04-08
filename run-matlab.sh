#!/bin/bash

#PBS -V
#PBS -N H3CRP
#PBS -l nodes=1:ppn=20,mem=60GB,walltime=48:00:00

module purge
module load intel/14.0.2
module load fftw/intel/3.3.4
module load matlab/2014a

cd /home/wang/matlab/h3-quantum-dynamics

j=(0 1 2 3) 
v=(0 1 2 3 4 5 6 7 8 9)

nj=${#j[@]};
nv=${#v[@]};

found=0
k=0
for j1 in ${j[@]}; do
    if [ $found -eq 1 ]; then break; fi
    for v1 in ${v[@]}; do
	if [ $found -eq 1 ]; then break; fi
	k=$((k+1))
	if [ $k -eq $PBS_ARRAYID ]; then
	    jRot=$j1
	    nVib=$v1
	    found=1
	fi
    done
done

matlab -nodesktop -nosplash -r "main($jRot, $nVib); exit" > H3-output-j${jRot}-v${nVib}.log 2>&1

exit



