% Investigate Galerkin discretisation of G

global N; N = 20; brackets

a0 = 2;
R = 3;
t = pi/2;
c0 = sum(nq*evan(-0.5*abs(a0).^2, a0), 2);
c0 = (exp(-1i*t*nhn).*c0)/norm(c0);
a = gample(c0, R);

plot(real(a), imag(a), '*k')
set(gca, 'DataAspectRatio', [1 1 1])

[ar, as] = meshgrid(a,a);
M = exp(-abs(ar-as).^2);
w = M \ (pi*ones(size(a)));

c = norm(exp(-0.5*abs(a(2))-0.5*abs(a)+conj(a)*a(2)));