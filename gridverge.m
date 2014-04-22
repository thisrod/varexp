% Fit cat states on a regular grid, and find the lagrange parameter at
% which the solution departs from an entire function, as the grid spacing
% varies.

global N; N = 20; brackets
c = sum(nq*evan([2i -2i], 'even'), 2);
hs = [0.3 1];
% hs = logspace(-1,0,3);
ls = logspace(-15,4,15);

nms = zeros(numel(hs), numel(ls));  rds = nms;
for i=1:numel(hs)
    figure
    [x,y] = meshgrid(-4:hs(i):4);  z = x(:)+1i*y(:);
    r = abs(z) <= 4;
    [M,f0] = evan(z(r), 'even');
    A = hs(i)*nq*M;
    for j=1:numel(ls)
        d = [A; ls(j)*eye(size(A))] \ [c; zeros(size(c))];
        nms(i,j) = norm(d);  rds(i,j) = norm(c-A*d);
        subplot(4,4,j)
        supplot([f0+log(d); z(r)])
    end
end

figure
for i=1:numel(hs)
    subplot(3,2,i);
    plot(rds(i,:), nms(i,:), '-k');
    title(num2str(hs(i)));
end