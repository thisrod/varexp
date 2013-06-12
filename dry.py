"""Dry run of the numerical experiment, integrating a quartic oscillator in a Fock basis by the explicit Euler formula."""

from numpy import arange, pi, zeros, abs, arange, concatenate, diag, exp, hstack, matrix, max, newaxis, sqrt, sum, vander
from scipy.misc import factorial
import matplotlib.pyplot as plt

N = 20
dt = 0.1
ts = dt*arange(100)
ats = zeros(ts.shape, dtype=complex)

nhn = matrix(diag(arange(N+1)*arange(-1,N)))
nq = matrix(diag(1./sqrt(factorial(arange(N+1)))))
aop = matrix(diag(sqrt(range(1,N+1)), 1))

def evan(f, a):
	expf = exp(f)[newaxis,:]
	va = vander(a, N+1)[:,::-1].T
	return matrix(expf*va)

a0 = 2
c = nq*evan([-0.5*2**2], [2])

for m in range(100):
	ats[m] = (c.H*aop*c)[0,0]
	c = c-1j*nhn*c*dt
	
plt.semilogy(ts, ats.real)
plt.show()

	