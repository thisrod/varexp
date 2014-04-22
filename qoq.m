% plot <a|psi(t)> for the quartic oscillator

global N; N = 30; n = (0:N)'; brackets
mesh = linspace(-5,5,50);  phasespace

H = n.^2-n;
c0 = nq*evan(2, 'even');
t = [0 0.04 0.1 0.16 0.18 0.2 0.24+0.1*rand() 1/2 2/3 0.9 0.96 1];
nt = numel(t);
f = 0:0.1:2*pi;

figure
set(gcf, 'DefaultAxesDataAspectRatio', [1 1 1])
for j = 1:nt
    c = exp(-0.5i*pi*t(j)*H).*c0;
    dc = -1i*H.*c;
    subplot(3,nt,j)
    zplot(mesh, mesh, Aps'*c);
    title(num2str(t(j)))
    subplot(3,nt,j+nt)
    zplot(mesh, mesh, Aps'*dc)
    subplot(3,nt,j+2*nt)
    zplot(mesh, mesh, DAps'*dc)
end