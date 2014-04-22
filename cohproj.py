"""This script compares a direct integration of the quartic oscillator by Euler's formula with one projected onto a superposition of coherent states at every step.
"""

import brackets; brackets.init(10); from brackets import *
from numpy import *
from numpy.linalg import lstsq, norm
from scipy.special import cbrt
import matplotlib.pyplot as plt

T = pi
c0 = nq*evan([-0.5*2**2], [2])
rcs = logspace(-4, 0, 20)
errs = zeros_like(rcs)
ndzs = zeros_like(rcs)

# set up initial guess
a0 = array([2, 1.5, 2-0.5j, 2.5, 2+0.5j])
f0 = -0.5*a0**2
c = sum(nq*evan(f0,a0), axis=1)
c = c/norm(c)
f0 = f0 - log(norm(c))
z0 = hstack((f0, a0))
dq = hstack((nq*evan(f0,a0), ndqr*evan(f0,a0)))

for i in range(rcs.size):
	dz = array(lstsq(dq, c0-c, rcond=rcs[i])[0]).flatten()
	ndzs[i] = norm(dz)
	z = z0 + dz
	d = sum(nq*evan(z[0:5],z[5:10]), axis=1)
	errs[i] = norm(d-c)
	
plt.figure(1)
plt.loglog(ndzs, errs)
plt.xlabel(r"$|\delta z|$")
plt.ylabel(r"$||\psi(z)\rangle-|\alpha\rangle|$")

plt.figure(2)
plt.loglog(rcs, errs)
plt.show()