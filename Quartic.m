% XSPDE version of variational quartic oscillator

% set up initial amplitudes

global N,  N = 10;  brackets
h = 0.6*sqrt(pi);  [Re Im] = ndgrid(h*(-10:10));
z = Re(:) + 1i*Im(:);  z = z(abs(z) < 3.9*h/sqrt(pi));
z = z + 2 + 0.1984 - 0.3435i;	% don't have a grid point at 2
R = numel(z);
[X, f0] = evan(z, 'even');
c = (nq*X)\(nq*evan(2,'even'));  f0 = f0 + log(c);
fprintf('initial residual = %.2e, component norm = %.2f\n', ...
	norm(nq*X*c - nq*evan(2,'even')), norm(c))
figure, quiver(real(z), imag(z), real(c), imag(c))
hold on, plot(2,0,'ok'), plot(z, '+k'), axis image
title 'initial variational weights'