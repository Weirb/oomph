import numpy as np
import matplotlib.pyplot as plt
from abc import abstractmethod
from enum import Enum


debug = False

def main():

	k = 6
	n = 2 ** k + 1
	h = 1 / (n + 1)

	a = 1
	b = a + 1

	x = np.linspace(a, b, n)

	poisson = False

	if poisson: # POISSON PROBLEM
		u_ex = x
		alpha = u_ex[0]
		beta = u_ex[-1]

		f = np.zeros(n)
		f[0] -= alpha / h ** 2
		f[-1] -= beta / h ** 2
			
		problem = MGPoissonProblem([a,b], f, k)

	else: # VARIABLE COEFFICIENT PROBLEM
		u0, u1 = 1, 0
		u_ex = u0*np.log(x) + u1

		alpha = u_ex[0]
		beta = u_ex[-1]

		f = np.zeros(n)
		f[0] -= alpha / h ** 2
		f[-1] -= beta * (1 / (b) / h + 1 / h ** 2)

		problem = MGVariableProblem([a,b], f, k)

	# Solve the problem
	solver = MGSolverEnum.TWO_GRID
	u, count = problem.solve(solver)
	
	# Log output
	residual = np.linalg.norm(problem.A(u)-f)
	abserror = np.linalg.norm(u_ex-u)
	relerror = abserror/np.linalg.norm(u_ex)

	print('Initial level : ', k)
	print('Solver        : ', solver)
	print('Residual      : ', residual)
	print('Abs error     : ', abserror)
	print('Rel error     : ', relerror)
	print('Iterations    : ', count)

	plt.plot(x, u)
	plt.plot(x, u_ex)
	plt.legend(["Multigrid", "Exact"])
	plt.show()



class MGBase:

	def __init__(self, f, k):
		self.f = f
		self.k = k
		self.n = 2**k+1
		self.u = np.zeros(self.n)
		self.jacobi_pre = 5
		self.jacobi_post = 5
		self.max_cycles = 300
	
	def prolong(self, v):
		n = len(v)
		return np.interp(np.linspace(0, 1, 2 * n - 1), np.linspace(0, 1, n), v)
	
	def restrict(self, v):
		return v[::2]
	
	def restrict_n(self, v, n):
		# Iteratively apply restriction to a vector to reach the desired level
		# Used in the mg F-cycle
		for i in range(n):
			v = self.restrict(v)
		return v
	
	def solve(self, solver_enum):

		if solver_enum == MGSolverEnum.TWO_GRID:
			solver = MGSolver.mg_two_grid
		elif solver_enum == MGSolverEnum.V_CYCLE:
			solver = MGSolver.mg_v_cycle
		elif solver_enum == MGSolverEnum.W_CYCLE:
			solver = MGSolver.mg_w_cycle
		elif solver_enum == MGSolverEnum.F_CYCLE:
			solver = MGSolver.mg_f_cycle
		else:
			print('Invalid solver selection.')
			return
		
		for count in range(self.max_cycles):
			u = solver(self)
			r = np.linalg.norm(self.A(u) - self.f)
			if r < 1e-5:
				break
		return u, count
	
	@abstractmethod
	def A(self, u):
		raise NotImplementedError
	
	@abstractmethod
	def create_mat(self, n):
		raise NotImplementedError
	
	@abstractmethod
	def smooth(self, u, f, max_iter):
		raise NotImplementedError


class MGSolverEnum(Enum):
	TWO_GRID = 1
	V_CYCLE = 2
	W_CYCLE = 3
	F_CYCLE = 4


class MGSolver(MGBase):

	def mg_two_grid(self):
		f = self.f
		# Pre-smooth
		u_fine = self.smooth(self.u, f, self.jacobi_pre)
		# Restrict the residual onto the coarser level
		residual = self.restrict(self.A(u_fine) - f)
		# Solve the problem exactly on the coarse level
		z = self.create_mat(int(np.ceil(self.n / 2)))
		u_coarse = np.linalg.solve(z, residual)
		# Apply the coarse grid correction
		u_fine = u_fine - self.prolong(u_coarse)
		# Post-smooth
		self.u = self.smooth(u_fine, f, self.jacobi_post)
		# Return final solution
		return self.u

	def mg_v_cycle(self):
		u = MGSolver._mg_cycle_runner(self, self.u, self.f, self.k, 1)
		self.u = u
		return u

	def mg_w_cycle(self):
		u = MGSolver._mg_cycle_runner(self, self.u, self.f, self.k, 2)
		self.u = u
		return u

	def _mg_cycle_runner(self, u, f, k, w):
		# w = 1 gives the v-cycle
		# w = 2 gives the w-cycle

		# Base case, solve the problem exactly
		if k == 3:
			if debug:
				print("solving on level {}".format(k))
			z = self.create_mat(2 ** k + 1)
			return np.linalg.solve(z, f)
		else:

			# Pre-smooth
			if debug:
				print("pre smoothing on level {}".format(k))
			u = self.smooth(u, f, self.jacobi_pre)

			# Compute the residual and restrict to coarser grid
			if debug:
				print("restriction on level {}".format(k))
			d = self.restrict(self.A(u) - f)
			v = np.zeros(np.shape(d))

			for s in range(1, w + 1):
				# Compute w recursive steps using the residual
				v = MGSolver._mg_cycle_runner(self, v, d, k - 1, w)

				# Add the correction back to the finer grid solution
			if debug:
				print("prolongation on level {}".format(k))
			u = u - self.prolong(v)

			# Post-smooth
			if debug:
				print("post smoothing on level {}".format(k))
			u = self.smooth(u, f, self.jacobi_post)

		if debug:
			print("returning value")
		return u

	def mg_f_cycle(self):
		n_max = int(np.log2(len(self.f)-1))
		u = MGSolver._mg_f_cycle_runner(self, self.f, 3, n_max)
		self.u = u
		return u

	def _mg_f_cycle_runner(self, f, n_min, n_max):
		# Perform a multigrid F-cycle, reaching the finest level n_max with
		# 2**n_max+1 elements and starting on the coarsest level n_min

		# Solve exactly on the coarsest level first
		z = self.create_mat(2 ** n_min + 1)
		u = np.linalg.solve(z, self.restrict_n(f, n_max - n_min))
		
		# Perform n V-cycles, moving up a level each iteration
		for i in range(n_min, n_max+1):

			# Interpolate to finer mesh first
			u = self.prolong(u)
			# Perform several V-cycles
			for k in range(20):
				u = MGSolver._mg_cycle_runner(self, u, self.restrict_n(f, n_max - i), i, 1)
				r = np.linalg.norm(self.A(u)-self.restrict_n(f,n_max-i))
				if r < 1e-8:
					break

		return u


class MGPoissonProblem(MGBase):

	def __init__(self, domain, f, k):
		super().__init__(f, k)
		self.a = domain[0]
		self.b = domain[1]

	def A(self, u):
		N = len(u)
		h = 1 / (N + 1)
		x = np.zeros(N)

		x[0] = 1 / h ** 2 * (-2 * u[0] + u[1])
		for i in range(1, N - 1):
			x[i] = 1 / h ** 2 * (u[i - 1] - 2 * u[i] + u[i + 1])
		x[-1] = 1 / h ** 2 * (u[-2] - 2 * u[-1])

		return x

	def create_mat(self, n):
		h = 1 / (n + 1)
		A = np.diag(-2 * np.ones(n))
		A += np.diag(np.ones(n - 1), 1)
		A += np.diag(np.ones(n - 1), -1)
		return A * 1 / h ** 2
	
	def smooth(self, u, f, max_iter):
		# Variables a and b are the domain end points
		# Weight parameter
		w = 2 / 3
		# Number of elements on the current level
		m = len(u)
		# Grid size on the current level
		h = 1 / (m + 1)
		# Array containing all of the interior points
		I = np.arange(1, m - 1)

		# Iterate until reached maximum number of iterations
		for k in range(max_iter):
			
			# Perform weighted Jacobi iteration
			# Treat end points separately
			# Vectorised operation to improve performance
			u[0] = (1 - w) * u[0] + w / 2 * (u[1] - h ** 2 * f[0])
			u[I] = (1 - w) * u[I] + w / 2 * (u[I - 1] + u[I + 1] - h ** 2 * f[I])
			u[-1] = (1 - w) * u[-1] + w / 2 * (u[-2] - h ** 2 * f[-1])
		return u


class MGVariableProblem(MGBase):
	def __init__(self, domain, f, k):
		super().__init__(f, k)
		self.a = domain[0]
		self.b = domain[1]
	
	def A(self, u):
		# Compute the matrix action for the variable coefficient discretisation

		# Number of elements on the current level
		N = len(u)
		# Grid size on the current level
		h = 1 / (N + 1)
		# Solution variable
		x = np.zeros(N)

		# Apply the discretisation stencil to the vector u
		# Treat end points separate from interior points
		x[0] = np.dot([-2 / h ** 2 - 1 / (self.a + h) / h, 1 / h ** 2 + 1 / (self.a + h) / h], u[:2])
		for i in range(1, N - 1):
			x[i] = np.dot(
				[
					1 / h ** 2,
					-2 / h ** 2 - 1 / (self.a + i * h) / h,
					1 / h ** 2 + 1 / (self.a + i * h) / h,
				],
				u[i - 1 : i + 2],
			)
		x[-1] = np.dot([1 / h ** 2, -2 / h ** 2 - 1 / (self.a + N * h) / h], u[-2:])

		return x

	def create_mat(self, N):
		# Create matrix of size NxN for solving the variable coefficient problem exactly

		# Grid size on level N
		h = 1 / (N + 1)
		# Matrix variable
		A = np.zeros((N, N))

		# Store the matrix entries into A
		A[0, :2] = [-2 / h ** 2 - 1 / (self.a + h) / h, 1 / h ** 2 + 1 / (self.a + h) / h]
		for j in range(1, N - 1):
			A[j, j - 1 : j + 2] = [
				1 / h ** 2,
				-2 / h ** 2 - 1 / (self.a + j * h) / h,
				1 / h ** 2 + 1 / (self.a + j * h) / h,
			]
		A[-1, N - 2 :] = [1 / h ** 2, -2 / h ** 2 - 1 / (self.a + N * h) / h]

		return A
	
	def smooth(self, u, f, max_iter):
		# Smoothing operation is done by the weighted Jacobi iteration

		# Variables a and b are the domain end points
		# Weight parameter
		w = 2 / 3
		# Number of elements on the current level
		m = len(u)
		# Grid size on the current level
		h = 1 / (m + 1)
		# Array containing all of the interior points
		I = np.arange(1, m - 1)

		# Iterate until reached maximum number of iterations
		for k in range(max_iter):
			# Perform weighted Jacobi iteration
			# Treat end points separately
			# Vectorised operation to improve performance
			u[0] = (1 - w) * u[1] - w / (-2 / h ** 2 - 1 / (self.a + h) / h) * (
				u[1] * (1 / h ** 2 + 1 / (self.a + h) / h) - f[0]
			)
			u[I] = (1 - w) * u[I] - w / (-2 / h ** 2 - 1 / (self.a + I * h) / h) * (
				u[I - 1] / h ** 2 + u[I + 1] * (1 / h ** 2 + 1 / (self.a + I * h) / h) - f[I]
			)
			u[-1] = (1 - w) * u[-1] - w / (-2 / h ** 2 - 1 / (self.a + m * h) / h) * (
				u[-2] / h ** 2 - f[-1]
			)

		return u




if __name__ == "__main__":
	main()
