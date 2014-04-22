% Plot the distance from a(|h>-|-h>) to |n>


N = 4;
[n,s,h] = ndgrid(0:N, logspace(-6,0,5), logspace(-4,0,4));

c = 2*exp(-0.5*h.^2).*h.^n./sqrt(factorial(n));
c(1:2:end,:,:) = 0;
c = reshape(c, [1 size(c)]);  c = repmat(c, N+1, 1);
d = reshape(eye(N+1), [N+1 N+1 1 1]);
d = repmat(d, 1, 1, numel(s), numel(h));


