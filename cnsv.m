% Sample a state by inverting G on random amplitudes, then look for
% unimportant singular directions of |Dpsi> along which the least singular
% value of G varies rapidly.

global N; N = 15; brackets
R = 5;
h = 0.01;
c0 = nq*evan(-2,2);

% Random sample amplitudes
a = 2+sqrt(0.5)*randn(R,2)*[1; 1i];
G = nq*evan(-0.5*abs(a).^2,a);
ef = pinv(G)*c0;
f = -0.5*abs(a).^2 + log(ef);

figure; subplot(1,2,1)
plot(real(a), imag(a), '*k');
set(gca, 'DataAspectRatio',[1 1 1])
for r = 1:R
    text(real(a(r)), imag(a(r)), sprintf('  %d',r))
end
title('Amplitudes')

[Ud,Sd,Vd] = svd([nq*evan(f,a) ndqr*evan(f,a)], 'econ');
[Ug,Sg,Vg] = svd(nq*evan(-0.5*abs(a).^2,a), 'econ');

ssdist = zeros(1,R);
for j = 1:R
    ssdist(j) = norm(Ud(:,1:j)*Ud(:,1:j)' - Ug(:,1:j)*Ug(:,1:j)');
end

figure; plot(1:R, ssdist, '+k')
xlabel('dimension of subspaces')
ylabel('operator distance between projectors')
title('Distances between ranges of low-rank approximations to G and D\psi')
    

ngs = zeros(1,2*R);     % Norms of grad sv
for p = 1:2*R
    da = h*Vd(R+1:2*R,p);
    Gp = nq*evan(-0.5*abs(a+da).^2, a+da);
    Gm = nq*evan(-0.5*abs(a-da).^2, a-da);
    A = [zeros(N+1), Gp-Gm; (Gp-Gm)', zeros(R)];
    ev = [Ug(:,end); Vg(:,end)];
    dsr = 0.5*ev'*A*ev;
    da = 1i*da;
    Gp = nq*evan(-0.5*abs(a+da).^2, a+da);
    Gm = nq*evan(-0.5*abs(a-da).^2, a-da);
    A = [zeros(N+1), Gp-Gm; (Gp-Gm)', zeros(R)];
    ev = [Ug(:,end); Vg(:,end)];
    dsi = 0.5*ev'*A*ev;
    ngs(p) = 0.5*norm([dsr, dsi])/h;
end

figure;
semilogy(1:2*R, diag(Sd), 'ok', 1:2*R, ngs, 'vk');
title('Derivatives of smallest sv of G along sds of D\psi')
legend('sv of D\psi', 'dsv of G')