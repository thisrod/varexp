% plots of quartic oscillator evolution
% idea: find the volume of phase space occupied by a ring of N oscillators with a Bose-Hubbard hamiltonian

% note spiral of zeroes at first local maximum.  Interference as larger amplitudes lap smaller ones.

% fucking Matlab pdf black-on-black labels


global N; N = 35; brackets
h = 0.1;  mesh = -5:h:5;  phasespace
ts = [0 0.05 1/7 1/6 2/7 1/5 1/2 1]/2;
ts = [ts 0.0342 0.0542 0.0585 0.0660 0.0763 0.0832 0.2002];
ts = sort(ts);

c0 = nq*evan(3.5, 'even');
for i = 1:numel(ts)
	t = ts(i);
	c = exp(2i*pi*nhn*t).*c0;	a = c'*aop*c;
	figure, zplot(mesh, mesh, Aps'*c), hold on, axis image
	plot(0, 0, 'ow', real(a), imag(a), '+w'), axis off
	saveTightFigure(sprintf('qrt%02d.pdf', i))

	figure, zplot(mesh, mesh, widi(c,z)), hold on, axis image
	plot(0, 0, 'ow', real(a), imag(a), '+w'), axis off
	saveTightFigure(sprintf('qrtw%02d.pdf', i))

	z0 = exp(-2i*pi*t*(2*abs(z).^2-1)).*z;
	figure, zplot(mesh, mesh, exp(-2*abs(z0-3.5).^2))
	axis image, axis off
	saveTightFigure(sprintf('qrtt%02d.pdf', i))
end

tt = 0:0.001:0.25;  m = numel(tt);
r = repmat(nan, 1, m);
w = repmat(nan, 1, m);

for j = 1:m
	c = exp(2i*pi*nhn*tt(j)).*c0;
%	o(j) = h^2*sum(abs(Aps'*c).^2)/pi;
	r(j) = abs(c'*aop*c);
	w(j) = h^2*sum(abs(Aps'*c).^2)^2/sum(abs(Aps'*c).^4);
end

ts = ts(ts <= 0.25);  wo = interp1(tt, w, ts);

figure, plot(tt, w/pi, '-k', ts, wo/pi, 'ok')
% title 'Coherent state spreading in quartic oscillator'
xlabel 't/\pi'
ylabel 'phase space area / \pi'
saveTightFigure qrta.pdf

figure, semilogy(tt, r, '-k')
title 'quadrature length'