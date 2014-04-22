% Examine the conditioning of expansion on regular and random grids,
% and the effect of parallel component conditioning

global N; N = 100; brackets
c = nq*evan(0.3+0.4i, 'even');
b = 2.5;      % radius of circular grid
R = 80;
nl = 15;
lz = logspace(-14,0,nl);
lw = logspace(-14,2,nl);

ax = -1.5*b:0.2:1.5*b;
[x,y] = meshgrid(ax);
zreg = x(:)+1i*y(:);

% set z to a regular grid
a = sqrt(2*pi/sqrt(3)/R)*b;
lim = ceil(2*b/sqrt(3)/a);
[i,j] = meshgrid(a*(-lim:lim));
z = i(:) + (0.5+0.5i*sqrt(3))*j(:);
z = z(abs(z)<=b);  nz = numel(z);
% set w to a random grid
w = gample(c, nz);

[Uz,Sz,~] = svd(nq*evan(z, 'even'), 'econ');
[Uw,Sw,~] = svd(nq*evan(w, 'even'), 'econ');
m = numel(diag(Sz));
    
figure(1);
subplot 121
zplot(ax, ax, (nq*evan(zreg, 'even'))'*c); hold on
plot(real(z), imag(z), '.w', 'MarkerSize', 12)
set(gca, 'DataAspectRatio', [1 1 1])
title([num2str(nz) ' points'])
subplot 122
zplot(ax, ax, (nq*evan(zreg, 'even'))'*c); hold on
plot(real(w), imag(w), '.w', 'MarkerSize', 12)
set(gca, 'DataAspectRatio', [1 1 1])

figure(2);
subplot 121
semilogy(1:m, diag(Sz), 'ok', 1:m, abs(Uz'*c), '+k', ...
    1:m, abs(Sz\Uz'*c), '*k')
title('Picard plot, regular grid'); xlabel n
legend('\sigma_n', 'u_n cpt', 'v_n cpt', 'Location', 'SouthWest')
subplot 122
semilogy(1:m, diag(Sw), 'ok', 1:m, abs(Uw'*c), '+k', ...
    1:m, abs(Sw\Uw'*c), '*k')
title('random grid'); xlabel n

Mz = nq*evan(z, 'even');  Mw = nq*evan(w, 'even');
Tz = zeros(nz*(nz-1), nz);
Tw = zeros(nz*(nz-1), nz);  k=1;
for i=1:nz
for j=1:i-1
    Tz(k,i) = Mz(:,j)'*c;
    Tz(k,j) = -Mz(:,i)'*c;
    Tw(k,i) = Mz(:,j)'*c;
    Tw(k,j) = -Mz(:,i)'*c;
    k = k+1;
end
end

Au = {Mz Mz Mw Mw};
Al = {eye(size(Mz)) Tz eye(size(Mw)) Tw};
l = [lz; lz; lw; lw];
nrhs = [N+1, nz*(nz-1), N+1, nz*(nz-1)];

nms = zeros(4, nl);  rds = nms;  cnp = nms; gs = nms;
for i = 1:4
for k=1:nl
    A = [Au{i}; l(i,k)*Al{i}];
    cnp(i,k) = cond(A);
    d = A \ [c; zeros(nrhs(i),1)];
    nms(i,k) = norm(d);  rds(i,k) = norm(c-M*d);
    p = abs(d)/sum(abs(d(:)));
    gs(i,k) = exp(-sum(p.*log(p)));
end
end

cz = Mz*Mz'*c;  cz = cz/norm(cz);  rz = norm(c-cz);
cw = Mw*Mw'*c;  cw = cw/norm(cw);  rw = norm(c-cw);

figure(3)
subplot 411
loglog(rds(1,:), nms(1,:), '-r', rds(2,:), nms(2,:), '-k', ...
    rds(3,:), nms(3,:), '--r', rds(4,:), nms(4,:), '--k')
xlabel('residual');  ylabel('solution norm')
subplot 412
loglog(l(1,:), rds(1,:), '-r', l(2,:), rds(2,:), '-k', ...
    l(3,:), rds(3,:), '--r', l(4,:), rds(4,:), '--k')
hold on
plot([l(1,1) l(1,nl)], [rz rz], ':k', ...
    [l(1,1) l(1,nl)], [rw rw], '-.k')
xlabel('lagrange parameter');  ylabel('residual')
subplot 413
loglog(l(1,:), cnp(1,:), '-r', l(2,:), cnp(2,:), '-k', ...
    l(3,:), cnp(3,:), '--r', l(4,:), cnp(4,:), '--k')
xlabel('lagrange parameter');  ylabel('condition of reg. problem')
subplot 414
semilogx(l(1,:), gs(1,:), '-r', l(2,:), gs(2,:), '-k', ...
    l(3,:), gs(3,:), '--r', l(4,:), gs(4,:), '--k')
xlabel('lagrange parameter');  ylabel('multiplicity of solution')

