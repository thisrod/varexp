% Evolve a coherent state for t=0.1, then fit a superposition of five
% coherent states to the result, truncated to an N=5 Fock basis.  Then
% expand the basis to N=10, and plot the singular values and the 
% components of vector from the approximate superposition to the state 
% along the singular directions in both psi and z space.  Finally, 
% plot distances from the shifted superposition to the target state, as
% more singular values are included in the shift, and as a smooth
% truncation parameter decreases.

global N
R = 5;
w = 1e-5;
t = 0.1;

% amplitude 2 coherent state and 5-component approximation
N = 7; brackets
c0 = exp(-1i*nhn*t).*(nq*evan(-2,2));

% 5-component approximation, randomly spread
a0 = 2+0.3*randn(1,R)+0.3i*randn(1,R);
f0 = -0.5*abs(a0).^2;
f0 = f0 + log(norm(c0)) - log(norm(sum(nq*evan(f0,a0),2)));
z = cohfit2(c0, [f0,a0], w);
f = z(1:R).'; a = z(R+1:2*R).';
abar = norm(a)/sqrt(R);


N = 10; brackets
c = exp(-1i*nhn*t).*(nq*evan(-2,2));
c0 = sum(nq*evan(f,a),2);
nwgt = diag(abs(sum(nq*evan(-0.5*abar^2,abar),2))); % Poisson
% nwgt = diag(exp(-(0:N)));     
% nwgt = eye(N+1);

figure(1)
plot(real(a), imag(a), '*k')
title(sprintf('residual = %.4f', norm(c-c0)))
set(gca, 'DataAspectRatio', [1 1 1])
xlabel('Re(\alpha)'); ylabel('Im(\alpha)')

figure(2)
[U,S,V] = svd(nwgt*[nq*evan(f,a), ndqr*evan(f,a)], 'econ');
dqcpt = U'*nwgt*(c-c0);
dzcpt = S\dqcpt;
dz = V*dzcpt;
semilogy(1:2*R, abs(dqcpt), '+b', 1:2*R, diag(S), 'ok', ...
    1:2*R, abs(dzcpt), '+r')
legend('\psi', 'S', 'z', 'Location', 'NorthWest')
title('Expansion of c-c_0')
xlabel('singular direction')
ylabel('component')

figure(3);
nms = zeros(1,R);
rsdls = zeros(1,R);
for lim=1:2*R
    Ut = U(:,1:lim); St = S(1:lim, 1:lim); Vt = V(:,1:lim);
    dz = Vt*(St\(Ut'*nwgt*(c-c0)));
    ff = f+dz(1:R).';  aa = a+dz(R+1:2*R).';
    cc = sum(nq*evan(ff,aa),2);
    nms(lim) = norm(cc);
    rsdls(lim) = norm(c-cc);
end
loglog(diag(S), rsdls, '-b', diag(S), nms, ':b'); hold on

slims = logspace(log10(min(diag(S))),log10(max(diag(S))),20);
nms = zeros(size(slims));
rsdls = zeros(size(slims));
for i=1:numel(slims)
    St = sqrt(S.^2+slims(i)^2*eye(2*R));
    dz = V*(St\(U'*nwgt*(c-c0)));
    ff = f+dz(1:R).';  aa = a+dz(R+1:2*R).';
    cc = sum(nq*evan(ff,aa),2);
    nms(i) = norm(cc);
    rsdls(i) = norm(c-cc);
end
loglog(slims, rsdls, '-g', slims, nms, ':g')
axis([min(slims), max(slims), 0.1, 1e3])
legend('trunc. error', 'trunc. norm', 'reg. error', 'reg. norm')
xlabel('bound on singular values')
title('Performance of truncation and regularisation')
