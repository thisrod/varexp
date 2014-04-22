% Follow the exact solution by nonlinear fitting
% now uses cohfit2 with w=1e-2, and includes number weighting

global N nwgt
N = 10; brackets
T = pi;
dts = 0.1;
R = 5;
w = 1e-5;

% amplitude 2 coherent state and 5-component approximation
c0 = nq*evan(-2,2);
% nwgt = diag(abs(c0).^2);
nwgt = eye(N+1);
nwgt(7:end,7:end) = 0;  % truncate

% 5-component approximation, randomly spread
a0 = 2+0.3*randn(1,R)+0.3i*randn(1,R);
f0 = -0.5*abs(a0).^2;
f0 = f0 + log(norm(c0)) - log(norm(sum(nq*evan(f0,a0),2)));

% exact solutions
rs = 0.5*pi*(1+nthroot(linspace(-1, 1, 100), 3));
crs = diag(c0)*exp(-1i*nhn*rs);
ars = diag(crs'*aop*crs);
figure(10)
plot(real(ars), imag(ars), '-k')

% propagation
for p=1:numel(dts)
    dt = dts(p);
    ts = (0:T/dt+1)*dt;
    as = zeros(size(ts)); as(:) = nan;
    rs = zeros(size(ts)); rs(:) = nan;
    cns = zeros(size(ts)); cns(:) = nan;


    f = f0; a = a0;
    for q=1:numel(ts)
        [z, rs(q)] = cohfit2(exp(-1i*nhn*ts(q)).*c0, [f,a], w);
        f = z(1:R).'; a = z(R+1:2*R).';
        ndq = [nq*evan(f,a), ndqr*evan(f,a)];
        c = sum(ndq(:,1:R),2);
        as(q) = c'*aop*c;
        if any(isnan(c)) || abs(as(q)) > 4
            break
        end
        cns(q) = cond(ndq);
        figure(3); subplot(5,7,q)
        plot(real(a), imag(a), '*k');
        figure(4); subplot(5,7,q)
        bar(sqrt(abs(sum([nq*evan(f,a),-exp(-1i*nhn*ts(q)).*c0],2)).^2));
    end
    figure(1)
    subplot(1,2,1)
    plot(real(ars), imag(ars), '-k', real(as), imag(as), '-b');
    set(gca, 'DataAspectRatio',[1 1 1])
    subplot(1,2,2)
    semilogy(ts, rs, '-k', ts, cns, ':k');
    
end