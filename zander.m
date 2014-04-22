function [A, p, q] = zander( z, n )
%zander Complex conjugate vandermonde
%   zander(z,n) where z is a complex number returns the row vector
%       [1 z z* z^2 zz* z^2* ...]
%   If z is a vector, it returns a matrix with rows of this form.
%   The extra return values are the powers of z* in p, and of z in q.
%   These are vectors.

z = z(:);  d = 0.5*(n+1)*(n+2);
A = zeros(numel(z), d);
p = zeros(1,d);  q = zeros(1,d);
p(1) = 0;  q(1) = 0;  A(:,1) = 1;  j = 1;
for i=1:n
    p(j+1:j+i+1) = 0:i;
    q(j+1:j+i+1) = i:-1:0;
    A(:,j+1:j+i) = diag(z)*A(:,j-i+1:j);
    j = j+i+1;
    A(:,j) = conj(z).^i;
end

