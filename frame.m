% Draw a normal ensemble of coherent states, and plot their frame
% projectors for a range of Tychonov regularisation parameters.

global N; N = 40; brackets
R = 15;
nl = 20;  tych = logspace(-5,1,nl);
pl = [20 18 15 13 10 1];
titles = {'components', 'dual ket', 'error', 'first three', '', '';
    '', 'weights', 'operator', 'dual kets', '', ''};

a = 2*sqrt(0.5)*[1 1i]*randn(2,R);
A = nq*evan(a,'even');

mesh = -5:0.2:5;
[x,y] = meshgrid(mesh);  z = x(:)+1i*y(:);
Aps = nq*evan(z,'even');

apops = zeros(nl, R, N+1);
enorms = zeros(1,nl);
wgts = zeros(nl, R);
conds = zeros(1, nl);
for i=1:nl
    M = pinv([A; tych(i)*eye(R)]);
    M = M(:,1:N+1);  apops(i,:,:) = M;
    wgts(i,:) = sqrt(abs(diag(M*M')));
    I = A*squeeze(apops(i,:,:));
    enorms(i) = sqrt(21) - norm(eye(21)-I(1:21,1:21), 'fro');
    conds(i) = cond(M);
end

% Labels

figure
for j=1:6
    subplot(6,7,7*j-6); axis off
    text(1, 0.5,titles{1,j}, ...
        'HorizontalAlignment', 'right', 'Units', 'Normalized')
    text(1, 0.3,titles{2,j}, ...
        'HorizontalAlignment', 'right', 'Units', 'Normalized')
end

% Plot ensembles and operators

for i=1:6
    l = pl(i);
    ty = squeeze(apops(l,:,:));
    
    subplot(6,7,i+1); hold on
    w = wgts(l,:)/max(wgts(l,:));
    for j=1:R
        plot(real(a(j)),imag(a(j)), 'ok', ...
            'MarkerFaceColor', hsv2rgb([0 0 1-w(j)]));
    end
    axis equal image;  if i>1; axis off; end
    
    subplot(6,7,7+i+1)    
    semilogy(1:R,wgts(l,:),'+k')
    
    subplot(6,7,2*7+i+1)    
    I = A*ty;
    zplot(0:15,0:15,eye(16)-I(1:16,1:16))
    axis image; if i>1; axis off; end
    
    dual = Aps'*squeeze(apops(l,:,:))';
    for j=1:3
        subplot(6,7,(2+j)*7+i+1)
        zplot(mesh, mesh, dual(:,j));
        axis off; hold on
        plot(real(a), imag(a), '.w')
        set(gca, 'DataAspectRatio', [1 1 1])
    end
end

% Plot singular vectors

figure
for i=1:6
	ty = squeeze(apops(pl(i),:,:));
    [U,S,V] = svd(ty, 'econ');
    for r=1:R
        subplot(6,R,R*i-R+r); hold on
        zplot(mesh, mesh, Aps'*V(:,r));
        axis image; hold on
        plot(real(a), imag(a), '.w')
        set(gca, 'DataAspectRatio', [1 1 1])
    end
end

figure
for i=1:6
	ty = squeeze(apops(pl(i),:,:));
    [U,S,V] = svd(ty, 'econ');
    subplot(5,7,i+1)
    semilogy(1:R,diag(S),'+k')
    
    subplot(5,7,7+i+1); hold on
    w = abs(U(:,1))/max(abs(U(:,1)));
    aw = 0.5*(1+angle(U(:,1))/pi);
    for j=1:R
        plot(real(a(j)),imag(a(j)), 'ok', ...
            'MarkerFaceColor', hsv2rgb([aw(j) 1 w(j)]));
    end
    axis equal image;  if i>1; axis off; end
        
    subplot(5,7,3*7+i+1); hold on
    w = abs(U(:,R))/max(abs(U(:,R)));
    aw = 0.5*(1+angle(U(:,R))/pi);
    for j=1:R
        plot(real(a(j)),imag(a(j)), 'ok', ...
            'MarkerFaceColor', hsv2rgb([aw(j) 1 w(j)]));
    end
    axis equal image;  if i>1; axis off; end
    
    subplot(5,7,2*7+i+1); hold on
    zplot(mesh, mesh, Aps'*V(:,1));
    axis image; hold on
    plot(real(a), imag(a), '.w')
    set(gca, 'DataAspectRatio', [1 1 1])

    subplot(5,7,4*7+i+1); hold on
    zplot(mesh, mesh, Aps'*V(:,R));
    axis image; hold on
    plot(real(a), imag(a), '.w')
    set(gca, 'DataAspectRatio', [1 1 1])
end

% L curve
figure
s = std(wgts.');
loglog(conds, enorms, '-k', conds(pl), enorms(pl), '.k', ...
    'MarkerSize', 12)
title('Performance of Tychonov conditioning')
xlabel('Condition of expansion operator')
ylabel('Quality of approximation')

