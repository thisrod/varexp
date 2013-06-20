N = 10; brackets
T = pi;
dts = logspace(-3,-1,6);

% amplitude 2 coherent state and 5-component approximation
c0 = nq*evan(-2,2);
a0 = 2+[0.5*1i.^(0:3), 0];  R = numel(a0);
f0 = -0.5*abs(a0).^2;
f0 = f0 + log(norm(c0)) - log(norm(sum(nq*evan(f0,a0),2)));
z = cohfit(c0, [f0,a0]); f0 = z(1:R); a0 = z(R+1:2*R);
c0 = sum(nq*evan(f0,a0),2);

% exact solution
rs = 0.5*pi*(1+nthroot(linspace(-1, 1, 100), 3));
crs = diag(c0)*exp(-1i*nhn*rs);
ars = diag(crs'*aop*crs);


h1 = figure;
set(gcf,'DefaultAxesDataAspectRatio',[1 1 1], ...
    'DefaultAxesVisible','off','DefaultAxesColor','none')
h2 = figure;

for p=1:numel(dts)
    figure(h1); subplot(2,3,7-p);
    dt = dts(p);
    ts = (0:T/dt)*dt;
    as = zeros(2, numel(ts)); as(:) = nan;
    rsdl = zeros(2, numel(ts)); rsdl(:) = nan;
    
    f = f0; a = a0; sty1 = 'g-';
    for q=1:numel(ts)
        c = sum(nq*evan(f,a),2);
        as(1,q) = c'*aop*c;
        dc = expm1(-1i*nhn*dt).*exp(-1i*nhn*ts(q)).*c0;
        ndq = [nq*evan(f,a), ndqr*evan(f,a)];
        if any(isnan(ndq(:))) || abs(as(1,q)) > 4
            as(1,q) = nan; sty1 = 'r-';
            plot(real(as(1,q-1)), imag(as(1,q-1)), 'or'); hold on
            break
        end
        dz = pinv(ndq, 0.1)*dc;
        f = f + dz(1:R).'; a = a + dz(R+1:2*R).';
        rsdl(1,q) = norm(dc-ndq*dz)/norm(dc);
    end
    
    f = f0; a = a0; sty2 = 'b-';
    for q=1:numel(ts)
        c = sum(nq*evan(f,a),2);
        as(2,q) = c'*aop*c;
        dc = expm1(-1i*nhn*dt).*c;
        ndq = [nq*evan(f,a), ndqr*evan(f,a)];
        if any(isnan(ndq(:))) || abs(as(2,q)) > 4
            as(2,q) = nan; sty2 = 'r-';
            plot(real(as(2,q-1)), imag(as(2,q-1)), 'or'); hold on
            break
        end
        dz = pinv(ndq, 0.1)*dc;
        f = f + dz(1:R).'; a = a + dz(R+1:2*R).';
        rsdl(2,q) = norm(dc-ndq*dz)/norm(dc);
    end
    
    plot(real(ars), imag(ars), '-k', real(as(1,:)), imag(as(1,:)), sty1, ...
        real(as(2,:)), imag(as(2,:)), sty2);
    if p==1; fmt = '\\tau=%.3f'; else fmt = '%.3f'; end
    text(0.65,1,sprintf(fmt, dt), ...
        'Units', 'Normalized', 'HorizontalAlignment', 'Center', ...
        'VerticalAlignment', 'Top')
    figure(h2); subplot(2,3,7-p);
    semilogy(ts,rsdl(1,:),sty1,ts,rsdl(2,:),sty2);
end
