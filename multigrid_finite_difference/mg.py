import numpy as np
import matplotlib.pyplot as plt

doc = False
poisson = False

a = 1
b = a + 1

if not poisson:
	A = lambda u: A_varcoef(u)
	create_mat = lambda n: create_mat_varcoef(n)
	smooth = lambda u, f, m: smooth_varcoef(u, f, m)
else:
	A = lambda u: A_poisson(u)
	create_mat = lambda n: create_mat_poisson(n)
	smooth = lambda u, f, m: smooth_poisson(u, f, m)


def main():

	max_cycles = 100

	
	k = 8
	n = 2 ** k + 1
	h = 1 / (n + 1)

	x = np.linspace(a, b, n)

	if not poisson:
		u_ex = np.log(x)

		alpha = u_ex[0]
		beta = u_ex[-1]

		f = np.zeros(n)
		f[0] -= alpha / h ** 2
		f[-1] -= beta * (1 / (a+n*h) / h + 1 / h ** 2)
	else:
		u_ex = x
		alpha = u_ex[0]
		beta = u_ex[-1]

		f = np.zeros(n)
		f[0] -= alpha / h ** 2
		f[-1] -= beta / h ** 2

	u = np.zeros(n)

	# Perform the f cycle
	# u = mg_f_cycle(f, 3, k)
	# r = np.linalg.norm(A(u) - f)
	# print(r)




	for i in range(max_cycles):
	    u = mg_two_grid(u, f)
	    # u = mg_cycle(u, f, k, 1)
	    # u = smooth(u, f, 2000)
	    r = np.linalg.norm(A(u) - f)
	    if r < 1e-5:
	        break

	u = np.linalg.solve(create_mat(n), f)
	print(np.linalg.norm(A(u)-f))
	print(np.linalg.norm(u_ex-u))

	# print(i, r)
	plt.plot(x, u)
	plt.plot(x, u_ex)
	plt.legend(["Multigrid", "Exact"])
	plt.show()


def mg_cycle(u, f, k, w):
	# w = 1 gives the v-cycle
	# w = 2 gives the w-cycle

	# Number of smoothing iterations
	jacobi_pre = 10
	jacobi_post = 10

	# Base case, solve the problem exactly
	if k == 3:
		if doc:
			print("solving on level {}".format(k))
		z = create_mat(2 ** k + 1)
		return np.linalg.solve(z, f)
	else:
		# Pre-smooth
		if doc:
			print("pre smoothing on level {}".format(k))
		u = smooth(u, f, jacobi_pre)

		# Compute the residual and restrict to coarser grid
		if doc:
			print("restriction on level {}".format(k))
		d = restrict(A(u) - f)
		v = np.zeros(np.shape(d))

		for s in range(1, w + 1):
			# Compute w recursive steps using the residual
			v = mg_cycle(v, d, k - 1, w)

			# Add the correction back to the finer grid solution
		if doc:
			print("prolongation on level {}".format(k))
		u = u - prolong(v)

		# Post-smooth
		if doc:
			print("post smoothing on level {}".format(k))
		u = smooth(u, f, jacobi_post)

	if doc:
		print("returning value")
	return u


def mg_f_cycle(f, n_max, n_min):
	# Perform a multigrid F-cycle, reaching the finest level n_max with
	# 2**n_max+1 elements and starting on the coarsest level n_min

	# Solve exactly on the coarsest level first
	z = create_mat(2 ** n_min + 1)
	u = np.linalg.solve(z, restrict_n(f, n_max - n_min))

	# Perform n V-cycles, moving up a level each iteration
	for i in range(n_min, n_max):
		# Interpolate to finer mesh first
		u = prolong(u)
		# Perform several V-cycles
		for k in range(20):
			u = mg_cycle(u, restrict_n(f, n_max - i), i, 1)
	return u


def mg_two_grid(u, f):
	jacobi_pre = 2
	jacobi_post = 2
	n = len(u)

	# Pre-smooth
	u_fine = smooth(u, f, jacobi_pre)
	# Restrict the residual onto the coarser level
	residual = restrict(A(u_fine) - f)
	# Solve the problem exactly on the coarse level
	z = create_mat(int(np.ceil(n / 2)))
	u_coarse = np.linalg.solve(z, residual)
	# Apply the coarse grid correction
	u_fine = u_fine - prolong(u_coarse)
	# Post-smooth
	u = smooth(u_fine, f, jacobi_post)
	# Return final solution
	return u


def smooth_poisson(u, f, max_iter):
	# Variables a and b are the domain end points
	# Weight parameter
	w = 2 / 3
	# Number of elements on the current level
	m = len(u)
	# Grid size on the current level
	h = 1 / (m + 1)
	# Array containing all of the interior points
	I = np.arange(1, m - 1)
	# Variable for tracking old solution
	# u_old = u

	# Iterate until reached maximum number of iterations
	for k in range(max_iter):
		# Perform weighted Jacobi iteration
		# Treat end points separately
		# Vectorised operation to improve performance
		u[0] = (1 - w) * u[0] + w / 2 * (u[1] - h ** 2 * f[0])
		u[I] = (1 - w) * u[I] + w / 2 * (u[I - 1] + u[I + 1] - h ** 2 * f[I])
		u[-1] = (1 - w) * u[-1] + w / 2 * (u[-2] - h ** 2 * f[-1])

		# # If we are no longer reducing the solution, stop iteration
		# if np.linalg.norm(u_old - u) < 1e-10:
		# 	break
		# 	# Keep the current solution for next time
		# u_old = u
	return u


def smooth_varcoef(u, f, max_iter):
	# Smoothing operation is done by the weighted Jacobi iteration

	# Variables a and b are the domain end points
	# Weight parameter
	w = 2 / 3
	# Number of elements on the current level
	m = len(u)
	# Grid size on the current level
	h = (b - a) / (m + 1)
	# Array containing all of the interior points
	I = np.arange(1, m - 1)
	# Variable for tracking old solution
	u_old = u

	# Iterate until reached maximum number of iterations
	for k in range(max_iter):
		# Perform weighted Jacobi iteration
		# Treat end points separately
		# Vectorised operation to improve performance
		u[0] = (1 - w) * u[1] - w / (-2 / h ** 2 - 1 / (a + h) / h) * (
			u[1] * (1 / h ** 2 + 1 / (a + h) / h) - f[0]
		)
		u[I] = (1 - w) * u[I] - w / (-2 / h ** 2 - 1 / (a + I * h) / h) * (
			u[I - 1] / h ** 2 + u[I + 1] * (1 / h ** 2 + 1 / (a + I * h) / h) - f[I]
		)
		u[-1] = (1 - w) * u[-1] - w / (-2 / h ** 2 - 1 / (a + m * h) / h) * (
			u[-2] / h ** 2 - f[-1]
		)

		# If we are no longer reducing the solution, stop iteration
		if np.linalg.norm(u_old - u) < 1e-10:
			break
			# Keep the current solution for next time
		u_old = u
	return u


def restrict(v):
	return v[::2]


def restrict_n(v, n):
	# Iteratively apply restriction to a vector to reach the desired level
	# Used in the mg F-cycle
	for i in range(n):
		v = restrict(v)
	return v


def prolong(v):
	n = len(v)
	return np.interp(np.linspace(0, 1, 2 * n - 1), np.linspace(0, 1, n), v)


def A_poisson(u):
	N = len(u)
	h = 1/ (N + 1)
	x = np.zeros(N)

	x[0] = 1 / h ** 2 * (-2 * u[0] + u[1])
	for i in range(1, N - 1):
		x[i] = 1 / h ** 2 * (u[i - 1] - 2 * u[i] + u[i + 1])
	x[-1] = 1 / h ** 2 * (u[-2] - 2 * u[-1])

	return x


def create_mat_poisson(n):
	h = 1 / (n + 1)
	A = np.diag(-2 * np.ones(n))
	A += np.diag(np.ones(n - 1), 1)
	A += np.diag(np.ones(n - 1), -1)
	return A * 1 / h ** 2


def A_varcoef(u):
	# Compute the matrix action for the variable coefficient discretisation

	# Number of elements on the current level
	N = len(u)
	# Grid size on the current level
	h = (b - a) / (N + 1)
	# Solution variable
	x = np.zeros(N)

	# Apply the discretisation stencil to the vector u
	# Treat end points separate from interior points
	x[0] = np.dot([-2 / h ** 2 - 1 / (a + h) / h, 1 / h ** 2 + 1 / (a + h) / h], u[:2])
	for i in range(1, N - 1):
		x[i] = np.dot(
			[
				1 / h ** 2,
				-2 / h ** 2 - 1 / (a + i * h) / h,
				1 / h ** 2 + 1 / (a + i * h) / h,
			],
			u[i - 1 : i + 2],
		)
	x[-1] = np.dot([1 / h ** 2, -2 / h ** 2 - 1 / (a + N * h) / h], u[-2:])

	return x


def create_mat_varcoef(N):
	# Create matrix of size NxN for solving the variable coefficient problem exactly

	# Grid size on level N
	h = (b - a) / (N + 1)
	# Matrix variable
	A = np.zeros((N, N))

	# Store the matrix entries into A
	A[0, :2] = [-2 / h ** 2 - 1 / (a + h) / h, 1 / h ** 2 + 1 / (a + h) / h]
	for j in range(1, N - 1):
		A[j, j - 1 : j + 2] = [
			1 / h ** 2,
			-2 / h ** 2 - 1 / (a + j * h) / h,
			1 / h ** 2 + 1 / (a + j * h) / h,
		]
	A[-1, N - 2 :] = [1 / h ** 2, -2 / h ** 2 - 1 / (a + N * h) / h]

	return A

if __name__ == "__main__":
	main()

