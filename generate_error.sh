#!/bin/bash

for k in {1..8}
do
	echo "Solving Helmholtz MG $k"
	main --linear_solver 1 --k_sq 100 --n_fourier 5 --min_ref $k \
	>> "logs/output_helmholtz_mg_$k.dat"

	echo "Solving Helmholtz SuperLU $k"
	main --linear_solver 0 --k_sq 100 --n_fourier 5 --min_ref $k \
	>> "logs/output_helmholtz_superlu_$k.dat"
done

