function d = discrepency(a,b,s)
	if nargin < 3, s = 2; end
	d = norm(a-b,s)/sqrt(norm(a,s)*norm(b,s));
end