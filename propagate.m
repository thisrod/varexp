function [ts, ms, rs] = propagate(dfun, mfun, z0, T, dt)
% Propagate the superposition z0 by applying dfun(t,c) to its Fock
% expansion, and projecting back to a superposition at each step.
% Return the results of applying mfun to the fock expansion, at
% each step.  Note that dfun is a differential, not a derivative.
% The results are passed to the continuations sfun or ffun as the
% integration succeeds or fails, and the result of that call returned.

global N nq ndqr

R = 0.5*numel(z0);
ts = (0:T/dt+1)*dt;
ms = zeros(numel(ts), numel(mfun(zeros(N+1,1))));
rs = zeros(numel(ts), 1);
f = z0(1:R).'; a = z0(R+1:2*R).';
for q=1:numel(ts)
    ndq = [nq*evan(f,a), ndqr*evan(f,a)];
    if any(isnan(ndq(:)))
        ts(q:end) = nan; ms(q:end) = nan; rs(q:end) = nan;
        return 
    end
    c = sum(ndq(:,1:R),2);
    ms(q,:) = mfun(c);
    dc = dfun(ts(q), c);
    dz = pinv(ndq, 0.1)*dc;
    f = f + dz(1:R).'; a = a + dz(R+1:2*R).';
    rs(q) = norm(dc-ndq*dz)/norm(dc);
end