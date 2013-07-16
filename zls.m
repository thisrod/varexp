function [ z ] = zls( X, b )
%zls Solve a least squares problem for a complex number and its conjugate
%   The problem solved is X*[z; conj(z)] = b

n = 0.5*numel(b);
Xl = X(:,1:n); Xr = X(:,n+1:end);
A = [real(Xl+Xr) real(1i*Xl-1i*Xr);...
    imag(Xl+Xr) imag(1i*Xl-1i*Xr)];
xy = A \ [real(b); imag(b)];
z = xy(1:n)+1i*xy(n+1:end);

end

