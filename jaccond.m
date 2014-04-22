% Fit a 5 component superposition to a coherent state, then investigate how
% the condition number of its derivative varies for small changes.

global N; N = 10; brackets
R = 5;
hs = linspace(-0.1,0.1,10);
w = 1.5e-3;
c0 = nq*evan(-2,2);

% Random sample amplitudes
a = 2+sqrt(0.5)*[1 1i]*randn(2,R);
G = nq*evan(-0.5*abs(a).^2,a);
ef = pinv(G)*c0;
f = -0.5*abs(a).^2 + log(ef).';

% amplitude 2 coherent state and 5-component approximation
disp(norm(sum([nq*evan(f,a) -c0], 2)))

figure; subplot(1,2,1)
plot(real(a), imag(a), '*k');
set(gca, 'DataAspectRatio',[1 1 1])
for r = 1:R
    text(real(a(r)), imag(a(r)), sprintf('  %d',r))
end
title('Amplitudes')

subplot(1,2,2)
bar(sqrt(sum(abs(nq*evan(f,a)).^2, 1)))
title('Weights')

basis = eye(R);
figure
for p = 1:R
    cns = zeros(2, numel(hs)); cns(:) = nan;
    for q = 1:numel(hs)
        b = a + hs(q)*basis(p,:);
        [U,S,V] = svd(nq*evan(-0.5*abs(b).^2,b), 'econ');
        cns(1,q)=S(end,end);
        b = a + 1i*hs(q)*basis(p,:);
        [U,S,V] = svd(nq*evan(-0.5*abs(b).^2,b), 'econ');
        cns(2,q)=S(end,end);
    end
    [U,S,V] = svd(nq*evan(-0.5*abs(a).^2,a), 'econ');
    b = a + 0.01*basis(p,:);  Gp = nq*evan(-0.5*abs(b).^2,b);
    b = a - 0.01*basis(p,:);  Gm = nq*evan(-0.5*abs(b).^2,b);
    A = [zeros(N+1), Gp-Gm; (Gp-Gm)', zeros(R)];
    ev = [U(:,end); V(:,end)];
    s = S(end,end);  ds = 0.5*ev'*A*ev;
    subplot(1,5,p)
    plot(hs, cns(1,:), '-k', [-0.1 0.1], s+[-0.1 0.1]*ds/0.02, ':k', ...
        hs, cns(2,:), '-r')
    if p==3
        title('Least singular value of perturbed expansion problem')
    end
end

figure
for p = 1:2*R
    cns = zeros(2, numel(hs)); cns(:) = nan;
    for q = 1:numel(hs)
        h = hs(q);
        z = [f a]; z(p) = z(p) + h;
        B = evan(z(1:R), z(R+1:2*R));
        cns(1,q)=cond([nq*B, ndqr*B]);
        z = [f a]; z(p) = z(p) + 1i*h;
        B = evan(z(1:R), z(R+1:2*R));
        cns(2,q)=cond([nq*B, ndqr*B]);
    end
    subplot(2,5,p)
    plot(hs, cns(1,:), '-k', hs, cns(2,:), '-r')
end