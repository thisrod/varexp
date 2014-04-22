% Fit 5 components to amplitude 2, with condition number weighting and a
% variety of weights.

global N
N = 10; brackets
ws = logspace(-3.5,0,20);

% amplitude 2 coherent state
c0 = nq*evan(-2,2);

% 5-component approximation, evenly spread
a0 = 2+[0.5*1i.^(0:3), 0]; R = numel(a0);
f0 = -0.5*abs(a0).^2;
f0 = f0 + log(norm(c0)) - log(norm(sum(nq*evan(f0,a0),2)));

cns = zeros(2,numel(ws));
rsdls = zeros(2,numel(ws));

figure(1)
for i = 1:numel(ws)
    z = cohfit3(c0, [f0,a0], ws(i));
    f = z(1:R).'; a = z(R+1:2*R).';
    l = nq*evan(f,a);
    rsdls(1,i) = norm(sum([l,-c0],2));
    ndq = [l, ndqr*evan(f,a)];
    cns(1,i) = cond(ndq);
    subplot(5,4,i)
    plot(real(a), imag(a), '*k')
end

% 5-component approximation, randomly spread
a0 = 2+0.3*randn(1,5)+0.3i*randn(1,5); R = numel(a0);
f0 = -0.5*abs(a0).^2;
f0 = f0 + log(norm(c0)) - log(norm(sum(nq*evan(f0,a0),2)));

figure(2)
for i = 1:numel(ws)
    z = cohfit3(c0, [f0,a0], ws(i));
    f = z(1:R).'; a = z(R+1:2*R).';
    l = nq*evan(f,a);
    rsdls(2,i) = norm(sum([l,-c0],2));
    ndq = [l, ndqr*evan(f,a)];
    cns(2,i) = cond(ndq);
    subplot(5,4,i)
    plot(real(a), imag(a), '*k')
end


figure(3); subplot(2,1,1)
loglog(cns(1,:), rsdls(1,:), '-k', cns(2,:), rsdls(2,:), '-g')
ylabel('residual')
xlabel('condition number')
subplot(2,1,2)
loglog(cns(1,:), ws, '-k', cns(2,:), ws, '-g')
xlabel('condition number')
ylabel('w')
