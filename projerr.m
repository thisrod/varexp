global N
N = 10; brackets
T = pi;
dts = logspace(-3,-1,6);
h1 = figure; h2 = figure;
set(h1,'DefaultAxesDataAspectRatio',[1 1 1], ...
    'DefaultAxesVisible','off','DefaultAxesColor','none')

% amplitude 2 coherent state and 5-component approximation
c0 = nq*evan(-2,2);
a0 = 2+[0.5*1i.^(0:3), 0];
f0 = -0.5*abs(a0).^2;
f0 = f0 + log(norm(c0)) - log(norm(sum(nq*evan(f0,a0),2)));
z0 = cohfit(c0, [f0,a0]);

% exact solutions
rs = 0.5*pi*(1+nthroot(linspace(-1, 1, 100), 3));
crs = diag(c0)*exp(-1i*nhn*rs);
ars = diag(crs'*aop*crs);

% computing
fmt = '%.3f';
for p=1:numel(dts)
    dt = dts(p);
    [t1, a1, r1] = propagate(...
        @(t,c)expm1(-1i*nhn*dt).*exp(-1i*nhn*t).*c0, ...
        @(c)c'*aop*c, ...
        z0, T, dt);
    a1(abs(a1)>4) = nan;
    [t2, a2, r2] = propagate(...
        @(t,c)expm1(-1i*nhn*dt).*c, ...
        @(c)c'*aop*c, ...
        z0, T, dt);
    a2(abs(a2)>4) = nan;
    figure(h1); subplot(2,3,7-p)
    if any(isnan(a1)); sty1='r-'; else sty1 = 'b-'; end
    if any(isnan(a2)); sty2='r-'; else sty2 = 'g-'; end
    plot(real(ars), imag(ars), '-k', real(a1), imag(a1), sty1 , ...
        real(a2), imag(a2), sty2);
    text(0.65,1,sprintf(fmt, dt), ...
        'Units', 'Normalized', 'HorizontalAlignment', 'Center', ...
        'VerticalAlignment', 'Top')
    figure(h2); subplot(2,3,7-p)
    semilogy(t1,r1,sty1,t2,r2,sty2)
end

