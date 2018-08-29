#!/bin/bash

kmax=6

# Get all of the multigrid data first
for ((k=1; k<=$kmax; k++))
do
	echo "Multigrid $k"

	h=`echo "scale=5;2^(1-$k)" | bc -l`
	file="logs/output_helmholtz_mg_$k.dat"

	# Get the error from each log
	# Only need to get this once for mg, not for superlu
	grep "Norm of error" $file \
	| awk -v z="$h" 'BEGIN{FS=":";ORS="\n"} {print z,$NF}' \
	>>"data/helmholtz_error.dat"

	# Get the time from the logs
	TIME=$(grep "Total time for linear solver" $file \
	| awk '{print $NF}')
	
	# Get the number of dofs from the logs
	DOFS=$(grep "Total time for linear solver" $file \
	| sed -r 's/[0-9]+\.[0-9]+//g' \
	| sed -r 's/[^0-9]*//g')
	
	echo $TIME $DOFS >>"data/helmholtz_mg_time.dat"

	# Get the number of iterations for FGMRES
	# This only applies to multigrid
	ITERATIONS=$(grep "Number of iterations to converge" $file \
	| sed -r 's/[^0-9]//g')

	echo $ITERATIONS $DOFS >>"data/helmholtz_mg_iter.dat"
done

# Now get the SuperLU data
for ((k=1; k<=$kmax; k++))
do
	echo "SuperLU $k"

	h=`echo "scale=5;2^(1-$k)" | bc -l`
	file="logs/output_helmholtz_superlu_$k.dat"

	# Get the time from the logs
	TIME=$(grep "Total time for linear solver" $file \
	| awk '{print $NF}')
	
	# Get the number of dofs from the logs
	DOFS=$(grep "Total time for linear solver" $file \
	| sed -r 's/[0-9]+\.[0-9]+//g' \
	| sed -r 's/[^0-9]*//g')
	
	echo $TIME $DOFS >>"data/helmholtz_superlu_time.dat"
done
