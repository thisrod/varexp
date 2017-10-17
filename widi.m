function W = widi(c, a)
%WIDI calculate the Wigner distribution for a pure quantum state
%
% WIDI(C, A) returns values of the Wigner distribution at the complex
% points A.  The vector C holds the number state components of the
% quantum state whose Wigner distribution is to be found.
%
% The algorithm is from pla-378-117.

W = zeros(size(a));

for m = 0:(length(c)-1)
	X = 2*(-1)^m/pi * exp(-2*abs(a).^2) .* polyval(LaguerreGen(m, 0), 4*abs(a).^2);
	W = W + abs(c(m+1))^2*X;
	for n = 0:(m-1)
		X = 2*(-1)^n/pi / prod(sqrt((n+1):m)) * exp(-2*abs(a).^2) ...
			.* (2*a).^(m-n) .* polyval(LaguerreGen(n, m-n), 4*abs(a).^2);
		W = W + 2*real(c(n+1)*conj(c(m+1))*X);
	end
end

end
