function B = norms(A,dim)
% like sum, but for the 2-norm

if nargin < 2, dim = 1; end
B = sqrt(sum(abs(A).^2,dim));
s = size(B);
B = reshape(B,[s(1:dim-1) s(dim+1:end)]);

end