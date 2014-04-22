% Investigate the stability of sinc reconstruction as the density reduces
% to the nyquist limit.


x = linspace(-5,5)';  k = 0:6;
xnq = (-5:pi/max(k):5)';
c = rand(numel(k),1);
f = sin(x*k)*c;  fnq = sin(xnq*k)*c;
f1 = sinc(0.5*max(k)*(x*ones(size(xnq'))-ones(size(x))*xnq'))*fnq;

figure
plot(x,f,'-k',x,f1,'--k'); hold on;
plot(xnq,fnq,'.r','MarkerSize',12)