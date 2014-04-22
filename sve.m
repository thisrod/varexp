% Investigate the singular values and vectors of the G operator.  Phase
% space is sampled on a grid with spacing 0.5, ket space is projected onto
% the fock states |0> through |10>.  The singular values are examined as
% the edges of the grid expand past alpha=sqrt(10).

global N; N = 24; brackets
mkrs = {'vk' 'sk' 'pk'};  mkr = 1;

h = 0.3;
for alim=4:6;
    msh = -alim:h:alim;
    [x,y] = meshgrid(msh);  x = x(:);  y = y(:);
    z = [-0.5*(x.^2+y.^2); x+1i*y];
    [U,S,V] = svd(h*sqrt(pi)*nq*evan(z)/pi);

    figure(1)
    semilogy(0:N, 1-diag(S), mkrs{mkr}, 0:N, ...
        1-sqrt(gammainc(1:N+1,alim^2)), '-k')
    mkr = mkr+1;
    hold on

    smpl = zeros(numel(msh));
    for i=1:N+1
        smpl(:) = V(:,i);
        figure(alim); subplot(5,5,i)
        zplot([-alim alim], [-alim alim], smpl)
        set(gca, 'DataAspectRatio', [1 1 1])
    end
end

figure(1)
set(gca, 'FontSize', 14)
title('Convergence of G as phase space boundary expands')
xlabel('n')
ylabel('1-\sigma_n')
legend('4x4', '6x6', '8x8', 'Location', 'SouthEast')
