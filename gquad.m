% Fit a superposition of coherent states to a cat, plot the Picard 
% and L curves, 
% interpolate phase space functions by gaussian polynomials, and find
% its norm.

global N; N = 40; brackets
mesh = -5:0.2:5;  phasespace


cat = sum(nq*evan([2i -2i], 'even'), 2);
R = 40;  a = gample(cat, R);
ls = logspace(-14,0,51);   % Lagrange parameters
deg = 4;    % degree of polynomial interpolants

[X, f0] = evan(a,'even');  A = nq*X;  [U,S,V]=svd(A, 'econ');
[Y, m, n] = zander(a, deg);  B = diag(exp(f0))*Y;
m = repmat(m,numel(m),1);  n = repmat(n,numel(n),1); 
G = pi*(m+n'==m'+n).*factorial(n+m');


ns = zeros(size(ls));  ms = zeros(size(ls));
rs = zeros(size(ls));  gs = zeros(size(ls));
ds = zeros(numel(zander(1,deg)), numel(ls));
fs = zeros(numel(a), 6);
for i=1:numel(ls)
    c = [A; ls(i)*eye(size(A))]\[cat; zeros(N+1,1)];
    ds(:,i) = B\c;  ms(i) = abs(ds(:,i)'*G*ds(:,i));
    ns(i) = norm(c);
    p = abs(c/ns(i)).^2;  gs(i) = exp(-sum(p.*log(p)));
    f = f0 + log(c);
    rs(i) = norm(cat - sum(nq*evan(f, a), 2));
    if any(i==(1:10:51)); fs(:,7-(i+9)/10) = f; end
end

figure
zplot(mesh, mesh, Aps'*cat);
axis image, hold on
for i=1:numel(a)
    plot(real(a(i)), imag(a(i)), 'ow')
end
saveTightFigure cample.pdf


figure
[x,y] = meshgrid(-2:0.3:2, -4:0.3:4);  z = x+1i*y;
[X, f1] = evan(z(:),'even');  bras = (nq*X)';
subplot 278
buf = zeros(size(z));  buf(:) = bras*cat;
zplot(-2:0.3:2, -4:0.3:4, buf);  set(gca, 'DataAspectRatio', [1 1 1])
for i = 1:6
	subplot(2,7,i+1); hold on
	supplot([fs(:,i); a])
	set(gca, 'DataAspectRatio', [1 1 1])
	title({sprintf('log(\\lambda) = %.1f', log10(ls(61-10*i))); ...
        sprintf('log(R) = %.1f', log10(rs(61-10*i)))})
    
    subplot(2,7,i+8)
	buf(:) = diag(exp(f1))*zander(z(:), deg)*ds(:,i);
    zplot(-2:0.3:2, -4:0.3:4, buf);  set(gca, 'DataAspectRatio', [1 1 1])
end

figure
nsv = size(diag(S));
semilogy(1:nsv, diag(S), 'ok', 1:nsv, abs(U'*cat), '+k', ...
    1:nsv, abs(S\U'*cat), '*k');
title(['Picard plot for expanding a cat over ' num2str(R) ...
    ' coherent states'])
legend('\sigma_r', '|<u_r|\psi>|', '|<u_r|\psi>|/\sigma_r', ...
    'Location', 'SouthWest')
xlabel 'r', ylabel 'singular value and vector sizes'
saveTightFigure capic.pdf

figure
imagesc(1:nsv, 0:N, abs(U).^2)
title('Components |<N|U>|^2 for singular kets of |A>')
set(gca, 'YDir', 'normal')
colormap gray; colormap(flipud(colormap))
xlabel('singular ket |U_r>');  ylabel('fock state |n>')

figure
title('L curve for Tychonov regularisation')
loglog(ns, rs, '-k')
xlabel('vector norm');  ylabel('residual')
title 'Cat state L curve'
saveTightFigure catyc.pdf

figure
subplot 312
loglog(ls, rs, '-k')
set(gca, 'XDir', 'Reverse')
xlabel('Lagrange multiplier');  ylabel('residual')

subplot 313
semilogx(ls, gs, '-k')
set(gca, 'XDir', 'Reverse')
xlabel('Lagrange multiplier');  ylabel('multiplicity')

figure
subplot 211
loglog(ls, ns, '-k')
set(gca, 'XDir', 'Reverse')
xlabel('Lagrange multiplier');  ylabel('vector norm')

subplot 212
loglog(ls, ms, '-k')
set(gca, 'XDir', 'Reverse')
xlabel('Lagrange multiplier');  ylabel('function norm')

