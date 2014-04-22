function as = gample( c, R )
%GAMPLE Sample amplitudes from the normalized entire distribution
%   cs is a vector of number state weights, R is the number of coherent
%   amplitudes to sample.  as is a column vector of coherent amplitudes.
%   This relies on a potentially dodgy assumption about how rapidy the
%   normalised entire expansion decays with amplitude.  This can be tested
%   with geck(cs).

global N nq

v = 5;  % assumes that |f(a)| <= exp(-|a|^2/v+0.5)
na = 5*R;    % Number of samples to take at once

as = zeros(R,1);  r = 0;
while r<R
    % [r R na]
    a = sqrt(v)*randn(na,2)*[1; 1i];
    p0 = exp(-0.5*abs(a).^2/v+0.5);
    ca = nq*evan(a, 'even');
    p = abs(ca'*c);
    a = a(rand(na,1) < p./p0);
    if numel(a) > R-r; a = a(1:R-r); end
    as(r+1:r+numel(a)) = a;
    r = r + numel(a);
end
