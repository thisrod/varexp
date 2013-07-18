% Sample a state by inverting G on random amplitudes, then look for
% unimportant singular directions of |Dpsi> along which the least singular
% value of G varies rapidly.

global N; N = 15; brackets
R = 5;
hs = -0.1:0.02:0.1;
c0 = nq*evan(-2,2);

% Random sample amplitudes
a = 2+sqrt(0.5)*randn(R,2)*[1; 1i];
G = nq*evan(-0.5*abs(a).^2,a);
ef = pinv(G)*c0;
f = -0.5*abs(a).^2 + log(ef);

figure
plot(real(a), imag(a), '*k');
set(gca, 'DataAspectRatio',[1 1 1])
for r = 1:R
    text(real(a(r)), imag(a(r)), sprintf('  %d',r))
end
title('Complex amplitudes')

[Ug,Sg,Vg] = svd(G, 'econ');
u = Ug(:,end); v = Vg(:,end); ev = [u; v];
sigma = Sg(end,end);

cns = zeros(2*R, numel(hs)); cns(:) = nan;
basis = [eye(R) 1i*eye(R)];
for p = 1:2*R
    for q = 1:numel(hs)
        da = hs(q)*basis(:,p);
        [U,S,V] = svd(nq*evan(-0.5*abs(a+da).^2, a+da), 'econ');
        cns(p,q)=S(end,end);
    end
end

gsym = zeros(1,2*R)
for p = 1:R
    x = real(a(p))*eye(N+1);  y = imag(a(p))*eye(N+1);
    gsym(p) = real(v(p)*u'*(aop'-x)*G(:,p));
    gsym(p+R) = real(v(p)*u'*(1i*aop'-y)*G(:,p));
end
   
figure
for p = 1:R
    subplot(2,R,p)
    plot(hs, cns(p,:), '-k', hs, cns(p+R,:), '-r', ...
        hs, sigma+gsym(p)*hs, ':k', hs, sigma+gsym(p+R)*hs, ':r')
    if p==1
        ylabel('\sigma_r')
    elseif p==3
        title('Variation of \sigma_r with x_r (black) and y_r (red)')
    end
    subplot(2,R,p+R)
    plot(hs, cns(p,:)-sigma-gsym(p)*hs, '-k', ...
        hs, cns(p+R,:)-sigma-gsym(p+R)*hs, '-r')
    if p==1; ylabel('\sigma_r-D\sigma_r dz'); end
end

