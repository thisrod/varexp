% plots of quartic oscillator evolution
% idea: find the volume of phase space occupied by a ring of N oscillators with a Bose-Hubbard hamiltonian

% note spiral of zeroes at first local maximum.  Interference as larger amplitudes lap smaller ones.

% fucking Matlab pdf black-on-black labels


global N; N = 35; brackets
h = 0.3;  mesh = -5:h:5;  phasespace

for a = [0 1 i 2+3i]
	c = aop'*nq*evan(a, 'even');
	figure, zplot(mesh, mesh, Aps'*c), hold on, axis image
end
%	saveTightFigure(sprintf('qrt%02d.pdf', i))
%	title(sprintf('Q.O. at t = %.2f*2\\pi', t))
