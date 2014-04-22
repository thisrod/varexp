% Solve the equation d|psi>/dt = |psi_0>-|psi>, where |psi_0> is a 
% quartic oscillator state, and |psi> is initially |psi_0> expanded over
% sampled amplitudes.  Monitor the convergence to |psi_0> and the condition
% of the Glauber operator.

global N; N = 15; brackets

t = pi/2;  c0 = exp(-1i*nhn*t).*(nq*evan(2, 'even'));
R = 15;  a = gample(c0, R);
le = 1e-9;  % Tychonov parameter for initial expansion
lv = 1e-4;  % variational problem
[A, f] = evan(a, 'even');  A = nq*A;
f = f + log([A'*A; le*eye(R)]\[A'*c0; zeros(R,1)]);
z0 = [f; a];

smax = 4;  ds = 0.03;  s = 0:ds:smax;  z = z0;
clear scapt zcapt
rls = zeros(size(s));
sings = zeros(size(s));
for i=1:numel(s)
    c = sum(nq*evan(z), 2);
    rls(i) = norm(c0-c);
    if i>2 && rls(i)>rls(i-1) && rls(i-2)>=rls(i-1)
        scapt = s(i); zcapt = z;
    end
    [U,S,V] = svd(evan(z(R+1:2*R), 'even'), 'econ');
    sings(i) = S(end:end);
    dq = [nq*evan(z) ndqr*evan(z)];
    zdot = [dq; lv*eye(size(dq))] \ [c0-c; 0*c];
    z = z+zdot*ds;
end

figure(1)
subplot 121
supplot(z0)
set(gca, 'DataAspectRatio', [1 1 1])
subplot 122
supplot(z)
set(gca, 'DataAspectRatio', [1 1 1])

figure(2)
subplot(2,1,1)
plot(s, rls, '-k')
ylabel('|c-c0|')
subplot(2,1,2)
plot(s, sings, '-k')
xlabel('s')
ylabel('lsv of |A>')