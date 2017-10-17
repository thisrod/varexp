function c = sqn(a,x)
%SQN expand squeezed states over number states
%
%SQN(A,X) where a is a coherent amplitude and x a squeezing factor returns a column vector of coefficients
%SQN(A,X) where a and x are n-vectors returns a matrix with n columns of expansion coefficients (NYI)
%
%See pra-13-2226, Equation 3.23

global N, n = (0:N)';

if (numel(a) > 1) || (numel(x) > 1)
	error 'vectors of squeezing parameters NYI'
end
if x==0
	error 'shortcut for zero NYI'
end

mu = cosh(x);  nu = sinh(x);
c = exp(-abs(a)^2/2 + conj(nu)*a^2/2/mu) * (nu/2/mu).^(n/2) ./ sqrt(mu*factorial(n));
for j = n'
	c(j+1) = c(j+1) * hermite(j, a/sqrt(2*mu*nu));
end

end