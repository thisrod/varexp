function X = evan(f,a)
% evan: Calculate the parts of Fock space brackets that depend on the
% amplitudes of the superposition.
% 
% X = evan(z)
% X = evan(f, a)
%
% npsi = sum(nq*evan(z),2)
% ndpsi = [nq*evan(z), ndqr*evan(z)]
%
% The amplitudes can be provided as a pair of vectors f and a, or a vector
% z = [f a].

global N

if nargin==1
    R = 0.5*numel(f);
    a = f(R+1:end); f = f(1:R);
end
a = a(:).';

X = ((ones(N+1,1)*a).^((0:N)'*ones(size(a))))*diag(exp(f));