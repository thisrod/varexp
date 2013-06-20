function [dc, dd, npsi, ndpsi] = check(f,a)

global N nq ndqr

npsi = zeros(N+1,1);
for m = 0:N
    for r = 1:numel(a)
        npsi(m+1) = npsi(m+1) + exp(f(r))*a(r)^m/sqrt(factorial(m));
    end
end
dc = norm(npsi-sum(nq*evan(f,a),2));

ndpsi = zeros(N+1,2*numel(a));
for m = 0:N
    for r = 1:numel(a)
        ndpsi(m+1,r) = exp(f(r))*a(r)^m/sqrt(factorial(m));
    end
end
for m = 1:N
    for r = 1:numel(a)
        ndpsi(m+1,r+numel(a)) = ...
            exp(f(r))*a(r)^(m-1)*sqrt(m/factorial(m-1));
    end
end
dd = norm(ndpsi-[nq*evan(f,a), ndqr*evan(f,a)]);