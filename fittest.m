% Experiment with fitting gaussian polynomials in z and z* to samples

D = 10;
x = -5:0.2:5; y = x;  [xs,ys] = meshgrid(x,y);  zs = xs(:)+1i*ys(:);
fs = sin(2*xs);

figure; subplot 241
V = diag(exp(-0.5*abs(zs).^2))*zander(zs, D);
for i = 0:D
    W = V(:,1:0.5*(i+1)*(i+2));
    c = W\fs(:);
    subplot(3,4,i+1); zplot(x, y, W*c);
    set(gca, 'DataAspectRatio', [1 1 1])
end
    