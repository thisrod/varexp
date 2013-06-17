"""Dry run of the numerical experiment, integrating a quartic
oscillator in a Fock basis by the explicit Euler formula."""

import brackets; brackets.init(10); from brackets import *
from numpy import *
from scipy.special import cbrt
import matplotlib.pyplot as plt

T = pi
c0 = array(nq*evan([-0.5*2**2], [2])).flatten()
dts = logspace(-4, 0, 20)
# dts = logspace(-3, 0, 3)

mrk = {2:"o", 4:"p", 7:"s", N:"*"}
sym = {}

rs = 0.5*T*(1+cbrt(linspace(-1, 1, 100)))
ct = exp(-1j*outer(nhn,rs))*c0[:,newaxis]
exact_ars = diag(matrix(ct).H*aop*matrix(ct))

splt = 0
for dt in dts:
	
	ts = dt*arange(ceil(T/dt))
	ars = zeros(ts.shape, dtype=complex)
	brs = zeros(ts.shape, dtype=complex)
	j = 0
	c = c0.copy()
	d = c0.copy()
	b = c0.copy()

	hsi = 1-1j*dt*nhn/(1+0.5*1j*dt*nhn)
	
	for m in range(ts.size):
		ars[m] = (matrix(d[:,newaxis]).H*aop*matrix(d[:,newaxis]))[0,0]
		brs[m] = (matrix(b[:,newaxis]).H*aop*matrix(b[:,newaxis]))[0,0]
		c = c-1j*nhn*c*dt	# Euler
		d = d/(1+1j*nhn*dt)	# backward Euler
		b = hsi*b			# semi-implicit
	
	plt.figure(1)
	for n in mrk.keys():
		plt.loglog(dt, abs(b[n]/c0[n]), "+b")
		plt.loglog(dt, abs(c[n]/c0[n]), mrk[n]+"r")
		sym[n], = plt.loglog(dt, abs(d[n]/c0[n]), mrk[n]+"k")
		
	plt.figure(2)
	plt.subplot(4,5,20-splt)
	plt.axis('off')
	plt.axis('equal')
	plt.plot(exact_ars.real, exact_ars.imag, 'k-', ars.real, ars.imag, 'r-', brs.real, brs.imag, 'b-')
	if splt==17:
		plt.title(r"$\tau=%.4f$" % dt)
	else:
		plt.title("$%.4f$" % dt)
	splt = splt + 1


plt.figure(1)
plt.title("Stiffness of Schrodinger's equation with $H=\hbar a^{2\dagger}a^2$")
plt.xlabel(r"Euler formula step $\tau$")
plt.ylabel("Growth of Fock coefficients $|c_n(4\pi)|/c_n(0)$")
plt.legend([sym[n] for n in sorted(mrk.keys())], ["$c_{%d}$" % n for n in sorted(mrk.keys())], loc="upper left")
plt.savefig("stiff.pdf")

plt.figure(2)
plt.savefig("becnvg.pdf")



	