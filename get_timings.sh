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

	FOLDER="logs/output_k_$k";
	FOLDER+="_n_$n";

	OUT_FOLDER="data/output_k_$k";
	OUT_FOLDER+="_n_$n";
	mkdir -p $OUT_FOLDER;

	# Get all of the multigrid data first
	for ((i=1; i<=$MAX_ITER; i++))
	do
		echo "Multigrid $i"
	
		h=`echo "scale=5;2^(1-$i)" | bc -l`
		file="$FOLDER/mg_$i.dat"
	
		# Get the error from each log
		# Only need to get this once for mg, not for superlu
		grep "Norm of error" $file \
		| awk -v z="$h" 'BEGIN{FS=":";ORS="\n"} {print z,$NF}' \
		>>"$OUT_FOLDER/error.dat"
	
		# Get the time from the logs
		TIME=$(grep "Total time for linear solver" $file \
		| awk '{print $NF}')
		
		# Get the number of dofs from the logs
		DOFS=$(grep "Total time for linear solver" $file \
		| sed -r 's/[0-9]+\.[0-9]+//g' \
		| sed -r 's/[^0-9]*//g')
		
		echo $TIME $DOFS >>"$OUT_FOLDER/mg_time.dat"
	
		# Get the number of iterations for FGMRES
		# This only applies to multigrid
		ITERATIONS=$(grep "Number of iterations to converge" $file \
		| sed -r 's/[^0-9]//g')
	
		echo $ITERATIONS $DOFS >>"$OUT_FOLDER/mg_iter.dat"
	done
	
	# Now get the SuperLU data
	for ((i=1; i<=$MAX_ITER; i++))
	do
		echo "SuperLU $i"
	
		h=`echo "scale=5;2^(1-$i)" | bc -l`
		file="$FOLDER/superlu_$i.dat"
	
		# Get the time from the logs
		TIME=$(grep "Total time for linear solver" $file \
		| awk '{print $NF}')
		
		# Get the number of dofs from the logs
		DOFS=$(grep "Total time for linear solver" $file \
		| sed -r 's/[0-9]+\.[0-9]+//g' \
		| sed -r 's/[^0-9]*//g')
		
		echo $TIME $DOFS >>"$OUT_FOLDER/superlu_time.dat"
	done
done
