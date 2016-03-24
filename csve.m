% Compute the the singular values and vectors of the G operator truncated
% to a disk, and compare to the calculated ones. 

global N; N = 24; brackets
mkrs = {'vk' 'sk' 'pk'};  mkr = 1;
l = zeros(size(mkrs));

h = 0.3;
for b=3:5;
    msh = -b:h:b;
    [x,y] = meshgrid(msh);  x = x(:);  y = y(:);
    a = x+1i*y;  in = abs(a)<=b; a = a(in);
    [U,S,V] = svd(h*pi^-0.5*nq*evan(a,'even'));

    figure(1)
    ll = plot(0:N, diag(S), mkrs{mkr}, ...
        0:N, sqrt(gammainc(b^2,(0:N)+1)), '-k');
    l(mkr) = ll(1);
    mkr = mkr+1;
    hold on

    smpl = zeros(numel(msh));
    for i=1:N+1
        smpl(in) = V(:,i);
        figure(b); subplot(5,5,i)
        zplot(msh, msh, smpl)
        set(gca, 'DataAspectRatio', [1 1 1])
    end

    figure(b), saveTightFigure(sprintf('rsv%d.pdf', b))
end

figure(1)
set(gca, 'FontSize', 14)
title('Convergence of |\Gamma_\chi> as phase space boundary expands')
xlabel('n')
ylabel('\sigma_n')
legend(l, 'b=3', '4', '5', 'Location', 'SouthWest')
saveTightFigure svgrid.pdf

figure(2)
set(gca, 'FontSize', 14)
v = 10.^-[6 4 2 1];  v = [v 0.1:0.2:0.9 1-v];
[n,m] = meshgrid(linspace(0, 24, 50));
[C,h] = contour(n,m,sqrt(gammainc(m,n+1)),v,'-k');
ylabel('b^2'); xlabel n
title('Singular values of |\Gamma_\chi>')
clabel(C,h)
saveTightFigure svcont.pdf