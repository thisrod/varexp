% Sample the vacuum state on regular and random sample amplitudes, and see
% how it converges.

global N; N = 30; brackets

c = nq*evan(0, 'even');
b = 2.5;  % radius of regular grid
nr = 10;  rs = logspace(1,3.5,nr);

% The truncated state the regular expansion should converge to

ct = c.*sqrt(gammainc(b^2, 1:N+1)).';

% Regular grid

rds = zeros(3, numel(rs));  hs = zeros(size(rs));
tsn = zeros(2, numel(rs));
for i=1:nr
    hs(i) = sqrt(pi/rs(i))*b;
    [x,y] = meshgrid(-b:hs(i):b);
    z = x(:)+1i*y(:);  z = z(abs(z) <= b);
    M = nq*evan(z, 'even');
    d = pi^-0.5*hs(i)*M'*c;
    T = M'*M*diag(d);   % kets are normalised
    tsn(1,i) = norm(T-T', 'fro');
    rds(1,i) = norm(c-pi^-0.5*hs(i)*M*d);
    rds(3,i) = norm(ct-pi^-0.5*hs(i)*M*d);
end

% Random grid, even weights, phases of entire expansion

tic; a = gample(c, round(rs(end))); toc
for i=1:nr
    r = round(rs(i));
    [M, f0] = evan(a(1:r), 'even');  M = nq*M;
    d = M'*c;  d = d/norm(d);
    f = f0 + log(d);  cfit = sum(nq*evan(f,a(1:r)), 2);
    T = M'*M*diag(d);
    tsn(2,i) = norm(T-T', 'fro');
    rds(2,i) = norm(c - cfit/norm(cfit));
end

figure
subplot 211
loglog(rs, rds(1,:), '-k', rs, rds(2,:), '--k', rs, rds(3,:), ':k')
xlabel('components');  ylabel('residual')
subplot 212
semilogx(rs, tsn(1,:), '-k', rs, tsn(2,:), '--k');
xlabel('components');  ylabel('tension')
print(gcf, 'converge.pdf', '-dpdf')
