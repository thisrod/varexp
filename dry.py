"""Dry run of the numerical experiment, integrating a quartic
oscillator in a Fock basis by the explicit Euler formula."""

import brackets; brackets.init(10); from brackets import *
from numpy import array, abs, arange, ceil, exp, logspace, pi, zeros
import matplotlib.pyplot as plt

T = 4*pi
a0 = 2
c0 = array(nq*evan([-0.5*2**2], [2])).flatten()

lbl = "$c_{%d}$"
mrk = {2:"ok", 4:"pb", 7:"sg", N:"*r"}

for dt in logspace(-4, 1, 20):

	ts = dt*arange(ceil(T/dt))
	c = c0.copy()
	for m in range(ts.size):
		c = c-1j*nhn*c*dt
	
	# FIXME: matplotlib mixes up legend unless orders match 
	for n in sorted(mrk.keys()):
		plt.loglog(dt, abs(c[n]/c0[n]), mrk[n], label=lbl % n)

plt.title("Stiffness of Schrodinger's equation with $H=\hbar a^{2\dagger}a^2$")
plt.xlabel(r"Euler formula step $\tau$")
plt.ylabel("Growth of Fock coefficients $c_n(4\pi)/c_n(0)$")
plt.legend([lbl % n for n in sorted(mrk.keys())], loc="upper left")
plt.savefig("dry.pdf")

	