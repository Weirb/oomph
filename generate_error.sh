#!/bin/bash

for k in {1..3}
do
	echo $k
	h=`echo "scale=5;2^(1-$k)" | bc -l`
	main --k_sq 10 --n_fourier 3 --min_ref $k \
	| grep "Norm of error" \
	| awk -v z="$h" 'BEGIN{FS=":";ORS="\n"} {print z,$NF}' \
	>>errors.dat
done
