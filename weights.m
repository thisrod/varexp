% Find weights for the best fit to exp(-x^2), as a sum of gaussians centred
% on x_r. 

R = 10;
x = linspace(-1,1,R).';
[rw, cl] = meshgrid(x);
A = exp(-0.5*(rw-cl).^2);
b = exp(-0.5*x.^2);
w0 = A \ ones(size(x));
v1 = A \ b;  w1 = v1./b;
AA = exp(-0.25*(rw.^2+cl.^2-6*rw.*cl));
bb = exp(-0.5*x.^2);
w2 = AA \ bb;
[x w0 w1 w2]