% Interpolate samples on a random grid

kmax = 5;   R = 20;
k = -kmax:kmax;
x = pi*(2*rand(1,R)-1).'; xg = (-pi:0.1:pi).';
A = exp(1i*x*k);    Ag = exp(1i*xg*k);
f = (x+pi).*(pi-x); fg = (xg+pi).*(pi-xg);
itp = Ag*pinv(A, 0.1);

figure
plot(x, zeros(size(x)), 'ok', x, f, '+k', xg, fg, ':k', xg, real(itp*f), '-k')

figure
for r=1:R
    subplot(4,5,r)
    plot(xg, real(itp(:,r)))
end