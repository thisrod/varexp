% Plot SVE of the frame operator |A>, where amplitudes are roots of 49.

global N; N = 100; brackets
mesh = linspace(-7,7,50);  phasespace
rs = 5:10:55;  nr = length(rs);

id = zeros(nr, N+1, N+1);
lsvs = zeros(nr, N+1, 5);
sws = zeros(nr,55);  sws(:) = nan;
nbar = sws;

for i=1:nr
    R = rs(i);
    h = 2*pi/R;  t = 4.5*exp((1:R)*1i*h);
    A = nq*evan(t, 'even');
    id(i,:,:) = A*pinv(A);
    [U,S,V] = svd(A,'econ');
    lsvs(i,:,:) = U(:,1:2*i-1:end);
    sws(i,1:R) = diag(S);
    nbar(i,1:R) = sum(diag(0:N)*abs(U).^2);
end

figure
for i=1:nr
    subplot(nr,1,i)
    plot(nbar(i,:),sws(i,:), 'ok')
    axis([0 55 0 3])
end
subplot(nr,1,nr), xlabel('n')
subplot(nr,1,3), ylabel('\sigma_n')

figure
set(gcf, 'DefaultAxesDataAspectRatio', [1 1 1])
for i=1:nr
    subplot(nr,1,i)
    zplot(0:50, 0:50, id(i,1:51,1:51))
end


figure
set(gcf, 'DefaultAxesDataAspectRatio', [1 1 1])
for i=1:nr
    for j = 1:5
        subplot(nr,6,6*i+j-5)
        zplot(mesh, mesh, Aps'*lsvs(i,:,j).')
    end
    subplot(nr,6,6*i-5), hold on
    plot(4.5*exp(2*pi*(1:rs(i))*1i/rs(i)), '.k')
    axis([-7 7 -7 7])
end
