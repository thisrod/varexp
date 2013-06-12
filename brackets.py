"""Calculate the required brackets for the quartic oscillator.

Throughout, the Fock basis comprises |0> through |N>.

The convention is that psi is a function of (psi0, psi1, ... psi_R, alpha0, alpha1, ..., alpha_R)

These definitions will generally be copied rather than imported, because of the difficulties of passing around N.
"""

from numpy import abs, arange, concatenate, diag, exp, hstack, matrix, max, newaxis, sqrt, sum, vander
from scipy.misc import factorial

def get_nhn(N): 
	"""Diagonal brackets of the quartic oscillator Hamiltonian between number states."""
	return arange(N+1)*arange(-1,N)

def get_nq(N):
	return matrix(diag(1./sqrt(factorial(arange(N+1)))))
def get_ndqr(N):
	return matrix(diag(sqrt(arange(1,N+1)/factorial(arange(N))), -1))
"""sum(nq*evan(f,a), axis=1) gives the bracket <N|psi>.
    hstack((nq*evan(f,a), ndqr*evan(f,a))) gives <N|Dpsi>
"""

def evan(f, a):
	expf = exp(f)[newaxis,:]
	va = vander(a, N+1)[:,::-1].T
	return matrix(expf*va)

	
def get_nhq(N):
	return matrix(diag(sqrt(concatenate(([0, 0], [m*(m-1)/factorial(m-2) for m in range(2,N+1)])))))
"""sum(nhq*evan(f,a), axis=1) gives the brackets of the Hamiltonian between number states and the superposition.
	"""

def get_aop(N):
	"Lowering operator"
	return matrix(diag(sqrt(range(1,N)), 1))

def check(f, a):
	"""Test routine.  Answers should be on the order of epsilon times the infinity norms of f and a."""
	R = len(f)
	x = hstack((nq*evan(f,a), ndqr*evan(f,a)))
	for m in range(N+1):
		for i in range(R):
			x[m,i] = x[m,i] - exp(f[i])*a[i]**m/sqrt(factorial(m))
	for m in range(1,N+1):
		for i in range(R):
			x[m,i+R] = x[m,i+R] - exp(f[i])*a[i]**(m-1)*sqrt(m/factorial(m-1))
	y = nhq*evan(f,a)
	for m in range(2,N+1):
		for i in range(R):
			y[m,i] = y[m,i] - exp(f[i])*a[i]**m*sqrt(m*(m-1)/factorial(m-2))
	return max(abs(x)), max(abs(y))