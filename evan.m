function X = evan(f,a)

global N

X = ((ones(N+1,1)*a).^((0:N)'*ones(size(a))))*diag(exp(f));