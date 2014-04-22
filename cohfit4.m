function [z, rsdl] = cohfit4(cs, z0, wgt)
% fit a superposition of coherent states to the fock expansion cs,
% starting from z.  A term w*log(cond(dq(z))) is added to the cost
% function in quadrature; there are no hard constraints.  Cohfit4 uses the
% global optimisation methods.

global nq ndqr

R = 0.5*numel(z0);
options = optimset('Display','off');
w0 = [real(z0), imag(z0)];
problem = createOptimProblem('fmincon', 'objective', @cost, 'x0', w0);
gs = GlobalSearch;
w = run(gs, problem);
z = (w(1:2*R)+1j*w(2*R+1:4*R)).';
f = z(1:R).'; a = z(R+1:2*R).';
rsdl = norm(sum([nq*evan(f,a),-cs],2));
                
                
function y = cost(w)
% Distance of z from state cs

f = w(1:R)+1j*w(2*R+1:3*R); a = w(R+1:2*R)+1j*w(3*R+1:4*R);
l = nq*evan(f,a);
ndq = [l, ndqr*evan(f,a)];
y = norm(sum([l,-cs],2))^2 + wgt^2*log(cond(ndq));

end

end