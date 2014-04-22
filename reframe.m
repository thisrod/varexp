function [f, a] = reframe(z)
% REFRAME: improve the condition of parameters describing a ket
%
% [f, a] = reframe(z)

global nq
z = z(:);  c = sum(nq*evan(z),2);
R = numel(z)/2;  f = z(1:R);  a = z(R+1:end);
figure;  plot(a, 'ok');  hold on;  axis equal

% Cull amplitudes with low weight
a = a(real(f) > -6);  R = numel(a);

% Add amplitudes between spread out points
as = repmat(a,1,R);  S = abs(as - as.') + diag(inf*ones(R,1));
[nndist, nn] = min(S);
j = nndist > 1 & nndist < 2;
a = [a; (a(j)+a(nn(j)))/2];  R = numel(a);

% Delete whichever of each pair has the lower index in a
as = repmat(a,1,R);  S = abs(as - as.') + diag(inf*ones(R,1));
[nndist, nn] = min(S);
a = a(nndist > 0.2 | nn < 1:R);  R = numel(a);

% Replace remaining isolated points with clusters
as = repmat(a,1,R);  S = abs(as - as.') + diag(inf*ones(R,1));
nndist = min(S);  j = nndist >= 2;
b = a(j);  h = rand(size(b));
a = [a(~j); b+0.3*exp(2i*pi*h); ...
    b+0.3*exp(2i*pi*(h+1/3)); b+0.3*exp(2i*pi*(h-1/3))];

% Resample in new frame
[M, f] = evan(a,'even');  A = nq*M;
f = f + log(pinv(A)*c);

plot(a, '+r')