% Evolve a coherent state for t=0.1, then fit a superposition of five
% coherent states to the result, truncated to an N=5 Fock basis.  Then,
% for the expansions N=6:10, plot the components of vector from the 
% superposition to the truncated state
% along the singular directions of |DQ>.

global N
R = 5;
w = 1e-5;
t = 0.1;

% amplitude 2 coherent state and 5-component approximation
N = 5; brackets
c = exp(-1i*nhn*t).*(nq*evan(-2,2));

% 5-component approximation, randomly spread
a0 = 2+0.3*randn(1,R)+0.3i*randn(1,R);
f0 = -0.5*abs(a0).^2;
f0 = f0 + log(norm(c0)) - log(norm(sum(nq*evan(f0,a0),2)));

figure(1)
plot(real(a), imag(a), '*k')
title(rsdl)

[z, rsdl] = cohfit2(c, [f0,a0], w);
f = z(1:R).'; a = z(R+1:2*R).';

for N=5:2:11
    brackets
    c = exp(-1i*nhn*t).*(nq*evan(-2,2));
    cprox = sum(nq*evan(f,a), 2);
    [U,S,V] = svd([nq*evan(f,a), ndqr*evan(f,a)]);
    incpt = U'*(c-cprox);
    outcpt = S\incpt;
    subplot(3,4,0.5*(N-3))
    semilogy(1:numel(incpt), abs(incpt), '+k')
    title(sprintf('N=%d', N))
    if N==5; ylabel('\psi components'); end
    subplot(3,4,0.5*(N-3)+4)
    semilogy(1:numel(outcpt), abs(outcpt), '+r')
    if N==5; ylabel('z components'); end
    subplot(3,4,0.5*(N-3)+8)
    semilogy(1:numel(diag(S)), diag(S), 'or')
    if N==5; ylabel('singular values'); end
end
    
