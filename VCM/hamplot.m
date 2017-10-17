% script to plot numerical features of a discretised Hamiltonian matrix
% assumes the following variables are set: ...

% Suffix convention: ev, ew, sw, u = lsv, v = rsv

Hew = nhn(1:R);
M = pinv(A)*diag(nhn)*A;
[ev,ew] = eig(M);  ew = diag(ew);
% ns = (norms(aop*A*ev)./norms(A*ev)).^2;
ns = sum(diag(0:N)*abs(A*ev).^2)./norms(A*ev).^2;
[ns, P] = sort(ns);
ev = ev(:,P);  ew = ew(P);
[Au,Asw,Av] = svd(A,'econ');  Asw = diag(Asw);
% Tychonov regularised H
Hy = pinv([A; laf*eye(R)])*[diag(nhn)*A; zeros(R)];
[Hyev, Hyew] = eig(Hy);  Hyew = diag(Hyew);
% yns = (norms(aop*A*Hyev)./norms(A*Hyev)).^2;
yns = sum(diag(0:N)*abs(A*Hyev).^2)./norms(A*Hyev).^2;
[yns, P] = sort(yns);
Hyev = Hyev(:,P);  Hyew = Hyew(P);

% Eigenvalues of discretised Hamiltonian and step operators

figure
zplot(1:R, 0:N, Au), axis image
title 'svs of A'

unit = eye(N+1,R);
Afock = pinv(A)*unit;
Afock = Afock/diag(norms(Afock));

figure, subplot 221
zplot(1:R, 1:R, ev), axis image
title 'evs of discretised H'

subplot 222
tmp2 = A*ev;  tmp2 = tmp2/diag(norms(tmp2));
imagesc(abs(tmp2(1:R,:)).^2), colormap gray, axis image
set(gca, 'YDir', 'normal')
title 'Fock discretised H'

subplot 223
zplot(1:R, 1:R, Afock), axis image
title 'number states'

subplot 224
tmp1 = A*Afock;  tmp1 = tmp1/diag(norms(tmp1));
imagesc(abs(tmp1(1:R,:)).^2), colormap gray, axis image
set(gca, 'YDir', 'normal')
title 'Fock number states'


figure
plot(0:R-1, yns, '.r', 0:R-1, ns, '.k', [0 R], [0 R], ':k')
title 'Average number in numerical eigenstates'
ylabel '<n>' , title 'excitation ns of discrete evs'
legend discrete regularised Location NorthWest

figure, subplot 311
plot(1:R, ew, 'vk', 1:R, Hew, '^k')
hold on, plot([1 R], 2/h*[1 1], '-k')
text(1, 2/h, 'Nyquist limit')
title 'H eigenvalues'
legend('discrete', 'exact', 'Location', 'SouthEast')

subplot 312
semilogy(1:R, Asw, '.k', [0 R], [laf laf]);
title 'Singular values of expansion'
xlabel n, ylabel '\sigma_n'

subplot 313
pew = Hew.*Asw.^2./(Asw.^2+laf.^2);

plot(1:R, Hyew, 'vk', 1:R, Hew, '^k', 1:R, pew, '.k')
hold on, plot([1 R], 2/h*[1 1], '-k')
title 'regularised H eigenvalues'
legend('discrete', 'exact', 'predicted', 'Location', 'NorthWest')
