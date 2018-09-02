#!/bin/bash

MAX_ITER=6

# Parameters we want the solution for
declare -a ks=(100 0 10 0);
declare -a ns=(5   0 0  5);

len=$(echo "${#ks[@]} - 1" | bc);

# Loop over the parameters we want to solve for
for j in `seq 0 $len`;
do
	k=${ks[$j]};
	n=${ns[$j]};

	echo "Solving problem k^2=$k, N=$n";

	FOLDER="logs/output_k_$k";
	FOLDER+="_n_$n";
	echo "Creating folder: $FOLDER";
	mkdir -p $FOLDER;

	# Compute the SuperLU solution
	for i in `seq 1 $MAX_ITER`;
	do
	
		echo "Solving SuperLU $i"
		main --linear_solver 0 --k_sq $k --n_fourier $n --min_ref $i \
		>> "$FOLDER/superlu_$i.dat"
	done
	
	# Compute the MG solution
	for i in `seq 1 $MAX_ITER`;
	do
		echo "Solving MG $i"
		main --linear_solver 1 --k_sq $k --n_fourier $n --alpha 0.5 --min_ref $i \
		>> "$FOLDER/mg_$i.dat"
	done
	echo;
done
