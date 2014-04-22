% Follow the exact solution by nonlinear fitting
% now uses cohfit5

% Parameters that seem to work:
% R=5; w=3e-5
% R = 20; w = 1e-4;     % convergence!

global N
N = 19; brackets
ts = 0.5*pi*[0:0.2:1];
R = 15;
w = 1e-4;
nsv = min(N+1,2*R);
hs = linspace(-0.1,0.1,5);


% amplitude 2 coherent state and 5-component approximation
c0 = nq*evan(-2,2);

% exact solutions
rs = 0.5*pi*(1+nthroot(linspace(-1, 1, 100), 3));
crs = diag(c0)*exp(-1i*nhn*rs);
ars = diag(crs'*aop*crs);

% propagation
as = zeros(size(ts)); as(:) = nan;
rs = zeros(size(ts)); rs(:) = nan;
cns = zeros(size(ts)); cns(:) = nan;

for q=1:numel(ts)
    [z, rs(q)] = cohfit5(exp(-1i*nhn*ts(q)).*c0, R, w);
    f = z(1:R).'; a = z(R+1:2*R).';
    ndq = [nq*evan(f,a), ndqr*evan(f,a)];
    c = sum(ndq(:,1:R),2);
    as(q) = c'*aop*c;
    cns(q) = cond(ndq);
    [U,S,V] = svd(ndq, 'econ');
    nbars = sum(((0:N).'*ones(1,nsv)).*abs(U).^2);
    
    figure(3); subplot(3,4,q)
    plot(real(a), imag(a), '*k')
    set(gca, 'DataAspectRatio',[1 1 1])
    figure(4); subplot(3,4,q)
    bar(sqrt(abs(sum([nq*evan(f,a),-exp(-1i*nhn*ts(q)).*c0],2)).^2));
    figure(5);  subplot(3,4,q)
    semilogy(1:nsv, abs(U'*(nhn.*c)), '+k', 1:nsv, diag(S), 'ok', ...
        1:nsv, abs(S\U'*(nhn.*c)), '*r')
    figure(6);  subplot(3,4,q)
    plot(1:nsv, nbars, 'xk')
    
    % plot changes of condition number
    figure(9+q)
    for p = 1:nsv
        cnhs = zeros(2, numel(hs)); cnhs(:) = nan;
        for q = 1:numel(hs)
            h = hs(q);
            zh = z + h*V(:,p);
            fh = zh(1:R).'; ah = zh(R+1:2*R).';
            cnhs(1,q) = cond([nq*evan(fh,ah), ndqr*evan(fh,ah)]);
            zh = z + 1i*h*V(:,p);
            fh = zh(1:R).'; ah = zh(R+1:2*R).';
            cnhs(2,q) = cond([nq*evan(fh,ah), ndqr*evan(fh,ah)]);
    end
    subplot(4,5,p)
    plot(hs, cnhs(1,:), '-k', hs, cnhs(2,:), '-r')
end
end

figure(1); subplot(1,2,1)
plot(real(ars), imag(ars), '-k', real(as), imag(as), '-b');
set(gca, 'DataAspectRatio',[1 1 1])
subplot(1,2,2)
semilogy(ts, rs, '-k', ts, cns, ':k');
