% Quartic oscillator with resampling

global N;  N = 30;  brackets

a = 2;  [~,f] = evan(a,'even');
dt = 1e-3;  T = 15*dt;
t = 0:dt:T;  S = numel(t);

[f, a] = reframe([f; a]);  z = [f; a];
c = sum(nq*evan(z),2);
ct = exp(1i*nhn*t).*repmat(c,1,S);

zt = cell(1,S);  cct = zeros(N+1,S);
for s=1:S
    zt{s} = z;  cct(:,s) = c;
    dz = 0;
    for j=1:4
        Dq = [nq*evan(z+dz/2), ndqr*evan(z+dz/2)];
        dz = pinv(Dq)*(1i*nhn.*c)*dt;
    end
    [f, a] = reframe(z + dz);  z = [f; a];
    c = sum(nq*evan(z),2);
end
zt{S} = z;  cct(:,S) = c;

% quadratures
Xt = diag(ct'*aop*ct);
XXt = diag(cct'*aop*cct);

figure
subplot 311
plot(Xt, ':k');  hold on;  plot(XXt, '-k');
set(gca, 'DataAspectRatio', [1 1 1])
title('Phase space trajectory')
subplot 312
plot(t, real(XXt), '-k', t, imag(XXt), '-r', ...
    t, real(Xt), ':k', t, imag(Xt), ':r')
title('Expected quadratures')
xlabel('t')
subplot 313
semilogy(t, diag(cct'*cct), '-k')
title('squared ket norm')


% components
figure
for s=1:S
    subplot(1,S,s)
    supplot(zt{s});  set(gca, 'DataAspectRatio', [1 1 1])
end