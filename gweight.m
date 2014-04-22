% Draw normally distributed points, and derive a set of weights such that
% the gaussian centred on each point can be integrated by a weighted sum of
% samples.

set(0, 'DefaultAxesFontSize', 14)

R = 15;  x = randn(R,1);  y = linspace(-2,2,R).';
[xr, xc] = meshgrid(x);  [yr, yc] = meshgrid(y);
Ax = exp(-(xr-xc).^2);  Ay = exp(-(yr-yc).^2);
rhs = ones(size(x))*sqrt(pi);
w = Ax\rhs;
v = Ay\rhs;

figure
subplot 211
plot(x, zeros(size(x)), 'ok', x, w, '+k')
title('Normally distributed abscissae')
xlabel('x_i')
ylabel('w_i')
subplot 212
plot(y, zeros(size(y)), 'ok', y, v, '+k')
title('Evenly spaced abscissae')
xlabel('x_i')
ylabel('w_i')

l = logspace(-10,1,20);      % reg. parameters
ry = zeros(size(l));  ny = zeros(size(l));
rx = zeros(size(l));  nx = zeros(size(l));
for i=1:numel(l)
    w = [Ax; l(i)*eye(R)] \ [rhs; zeros(size(x))];
    nx(i) = norm(w);  rx(i) = norm(Ax*w-rhs);
    v = [Ay; l(i)*eye(R)] \ [rhs; zeros(size(y))];
    ny(i) = norm(v);  ry(i) = norm(Ay*v-rhs);
end

figure
subplot 311
loglog(nx, rx, ':k', ny, ry, '-k');
legend('random', 'regular', 'Location', 'NorthEast')
xlabel('solution norm'); ylabel('residual')
title('Tychonov curves')
subplot 312
loglog(l, rx, ':k', l, ry, '-k');
xlabel('parameter'); ylabel('residual')
subplot 313
loglog(l, nx, ':k', l, ny, '-k');
xlabel('parameter'); ylabel('solution norm')
