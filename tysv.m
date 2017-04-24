s = [logspace(-3,-0.5,5) logspace(-0.5,0.5,10) logspace(0.5,3,5)];
a = exp(log(10)*(2*rand(30,1)-3));

figure
set(gcf, 'DefaultLineMarkerSize', 12, ...
    'DefaultAxesPlotBoxAspectRatio', [1.1 1 1], ...
    'DefaultAxesXTick', 10.^[-3 0 3])
for i=1:3
    subplot(1,3,i)
    loglog(s, s./(s.^2+1), '-k', a, a./(a.^2+1), '.k')
    axis([1e-3 1e3 1e-3 1])
    if i==2; xlabel('\sigma/\lambda'); end
    if i==1
        title('Too regularised')
        ylabel('\lambda\sigma^+')
    end
    if i==3 title('Under regularised'); end
    if i>1; set(gca, 'YTick', []); end
    a = 100*a;
end

saveTightFigure tysv.pdf