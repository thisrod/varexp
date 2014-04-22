% Expand the exact solutions over a grid, using the generating function.
% Find Picard plots for expanding their derivatives over the same grid,
% and for fitting the grid variationally to the derivatives.

global N; N = 25; brackets

set(0, 'DefaultAxesFontSize', 14)

h = 0.3;  alim = 6;
c0 = nq*evan(-2,2);
ts = linspace(0, 0.5*pi, 12);
[x,y] = meshgrid(-alim:h:alim);
a = x(:)+1i*y(:);
zs = [-0.5*abs(a).^2; a];
xpn = h*evan(zs)'*diag(factorial(0:N).^-0.5); % backward map

figure(1)
semilogy(0:N, abs(c0), 'vk', 0:N, abs(nhn.*c0), 'pk')
xlabel('Particle number n')
ylabel('|<n|?>|')
legend('|\psi>', 'H|\psi>')
title('Expansions over number states')

smpl = zeros(numel(-alim:h:alim));
errs = zeros(2,numel(ts));
for j=1:numel(ts)
    c = exp(-1j*nhn*ts(j)).*c0;  Hc = nhn.*c;
    gs = xpn*(exp(-1j*nhn*ts(j)).*c0);
    f = -0.5*abs(a).^2 + log(h*gs/pi);
    ndq = [nq*evan(f,a), ndqr*evan(f,a)];  [U,S,V] = svd(ndq, 'econ');
    
    figure(2); subplot(4,3,j)
    smpl(:) = gs;  zplot([-alim alim], [-alim, alim], smpl);
    set(gca, 'DataAspectRatio', [1 1 1])
    if j==2; title('State from time 0 to \pi/2'); end
    
    figure(3); subplot(4,3,j)
    smpl(:) = xpn*Hc;  zplot([-alim alim], [-alim, alim], smpl);
    set(gca, 'DataAspectRatio', [1 1 1])
    if j==2; title('Derivative from time 0 to \pi/2'); end
    
    figure(4); subplot(4,3,j)
    smpl(:) = xpn*U(:,1);  zplot([-alim alim], [-alim, alim], smpl);
    set(gca, 'DataAspectRatio', [1 1 1])
    if j==2; title('Largest singular ket of variational problem'); end   
    
    figure(5); subplot(4,3,j)
    smpl(:) = xpn*U(:,end);  zplot([-alim alim], [-alim, alim], smpl);
    set(gca, 'DataAspectRatio', [1 1 1])
    if j==2; title('Smallest singular ket of variational problem'); end
    
    figure(6); subplot(4,3,j)
    vHc = ndq*V*(S\U'*Hc);
    errs(1,j) = norm(vHc-Hc);
    smpl(:) = xpn*vHc;  zplot([-alim alim], [-alim, alim], smpl);
    set(gca, 'DataAspectRatio', [1 1 1])
    if j==2; title('Variational derivative'); end
    
    figure(8); subplot(4,3,j)
    semilogy(0:N, diag(S), 'ok', 0:N, abs(U'*Hc), '+k', ...
        0:N, abs(S\U'*Hc), 'pk')
    if j==2
        title('Picard plots for variational problem')
    end
end

figure(9)
plot(ts, errs(1,:), 'ok')
xlabel('time')
ylabel('variational error')
