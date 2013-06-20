"""Initialise the required brackets for the quartic oscillator.

Throughout, the Fock basis comprises |0> through |N>.

The convention is that psi is a function of (psi0, psi1, ... psi_R, alpha0, alpha1, ..., alpha_R)

These definitions will generally be copied rather than imported, because of the difficulties of passing around N.
"""


import numpy as _np, scipy.misc as _sp

def evan(f, a):
	expf = _np.exp(f)[_np.newaxis,:]
	va = _np.vander(a, N+1)[:,::-1].T
	return _np.matrix(expf*va)


def init(N):
	global nhn, nq, ndqr, nhq, aop

	globals()["N"] = N

	"""Diagonal brackets of the quartic oscillator Hamiltonian between number states."""
	nhn = _np.arange(N+1)*_np.arange(-1,N)

	"""sum(nq*evan(f,a), axis=1) gives the bracket <N|psi>.
	    hstack((nq*evan(f,a), ndqr*evan(f,a))) gives <N|Dpsi>
	"""
	nq = _np.matrix(_np.diag(1./_np.sqrt(_sp.factorial(_np.arange(N+1)))))
	ndqr = _np.matrix(_np.diag(_np.sqrt(_np.arange(1,N+1)/_sp.factorial(_np.arange(N))), -1))

	"""sum(nhq*evan(f,a), axis=1) gives the brackets of the Hamiltonian between number states and the superposition.
	"""
	nhq = _np.matrix(_np.diag(_np.sqrt(_np.concatenate(([0, 0], [m*(m-1)/_sp.factorial(m-2) for m in range(2,N+1)])))))

	"Lowering operator"
	aop = _np.matrix(_np.diag(_np.sqrt(_np.arange(1,N+1)), 1))


def _check(f, a):
	"""Test routine.  Answers should be on the order of epsilon times the infinity norms of f and a."""
	R = len(f)
	x = _np.hstack((nq*evan(f,a), ndqr*evan(f,a)))
	for m in range(N+1):
		for i in range(R):
			x[m,i] = x[m,i] - _np.exp(f[i])*a[i]**m/_np.sqrt(_sp.factorial(m))
	for m in range(1,N+1):
		for i in range(R):
			x[m,i+R] = x[m,i+R] - _np.exp(f[i])*a[i]**(m-1)*_np.sqrt(m/_sp.factorial(m-1))
	y = nhq*evan(f,a)
	for m in range(2,N+1):
		for i in range(R):
			y[m,i] = y[m,i] - _np.exp(f[i])*a[i]**m*_np.sqrt(m*(m-1)/_sp.factorial(m-2))
	return _np.max(_np.abs(x)), _np.max(_np.abs(y))