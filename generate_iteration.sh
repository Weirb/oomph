#!/bin/bash

MAX_ITER=2

# Parameters we want the solution for
declare -a ks=(0 10 20 50 100);
declare -a ns=(0 5);

len=$(echo "${#ks[@]} - 1" | bc);

# Loop over the n parameters
for p in `seq 0 1`;
do
	# Loop over the k paramters
	for j in `seq 0 $len`;
	do
		k=${ks[$j]};
		n=${ns[$p]};
	
		echo "Solving problem k^2=$k, N=$n";
	
		FOLDER="logs/iteration_k_$k";
		FOLDER+="_n_$n";
		echo "Creating folder: $FOLDER";
		mkdir -p $FOLDER;
			
		# Compute the MG solution
		for i in `seq 1 $MAX_ITER`;
		do
			echo "Solving MG $i"
			main --linear_solver 1 --k_sq $k --n_fourier $n --alpha 0.5 --min_ref $i \
			>> "$FOLDER/mg_$i.dat"
		done
		echo;
	done
done
