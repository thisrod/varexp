function [z, rsdl] = cohfit5(cs, R, wgt)
% fit a superposition of coherent states to the fock expansion cs,
% starting from a sampling of the glauber analytic expansion.
% A term w*log(cond(dq(z))) is added to the cost
% function in quadrature; there are no hard constraints.

global nq ndqr

N = numel(cs)-1;  Nbar = (0:N)*abs(cs).^2;

% initial guess by rejection sampling
a = [];  f = [];
while numel(a)<R
    amp = sqrt(Nbar)*[4 4i]*(rand(2,1)-0.5);
    gmp = exp(-0.5*abs(amp)^2)*conj(amp).^(0:N)./sqrt(factorial(0:N))*cs;
    if rand()>abs(gmp); continue; end
    a = [a amp];  f = [f -0.5*abs(amp)^2+1i*angle(gmp)];
end
z0 = [f a];

options = optimset('Display','notify');
w0 = [real(z0), imag(z0)];
w = fminunc(@(w)cost(w),w0,options);
z = (w(1:2*R)+1j*w(2*R+1:4*R)).';
f = z(1:R).'; a = z(R+1:2*R).';
rsdl = norm(sum([nq*evan(f,a),-cs],2));
                
                
function y = cost(w)
% number weighted distance of z from state cs

f = w(1:R)+1j*w(2*R+1:3*R); a = w(R+1:2*R)+1j*w(3*R+1:4*R);
l = nq*evan(f,a);
ndq = [l, ndqr*evan(f,a)];
y = norm(sum([l,-cs],2))^2 + wgt^2*log(cond(ndq));

end

end