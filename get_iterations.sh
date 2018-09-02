#!/bin/bash

MAX_ITER=2;

# Parameters we want the solution for
declare -a ks=(0 10 20 50 100);
declare -a ns=(0 5);

len=$(echo "${#ks[@]} - 1" | bc);

for i in `seq 1 $MAX_ITER`;
do
	# Loop over the n
	for p in `seq 0 1`;
	do
		n=${ns[$p]};
		OUT_FILE="data/iteration_n_$n";
	
		# Loop over the k
		for j in `seq 0 $len`;
		do
			k=${ks[$j]};
		
			FOLDER="logs/iteration_k_$k";
			FOLDER+="_n_$n";
		
			echo "Multigrid $i k=$k"
			
			file="$FOLDER/mg_$i.dat"
			
			# Get the number of dofs from the logs
			DOFS=$(grep "Total time for linear solver" $file \
			| sed -r 's/[0-9]+\.?[0-9]*$//g' \
			| sed -r 's/[^0-9]*//g')
			
			# Get the number of iterations for FGMRES
			ITERATIONS=$(grep "Number of iterations to converge" $file \
			| sed -r 's/[^0-9]//g')
			
			echo -n $ITERATIONS, >>"$OUT_FILE";
		done
		echo $DOFS>>$OUT_FILE;
	done
done
