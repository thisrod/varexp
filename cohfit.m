function [z, rsdl] = cohfit(cs, z0, ~)
% fit a superposition of coherent states to the fock expansion cs,
% starting from z.

options = optimset('Display','off','Algorithm','active-set');
[z, rsdl] = fmincon(@(z)cost(z),z0,...
                    [],[],[],[],[],[],@constraint,options);
z = z.';
                
                
function y = cost(z)
% Distance of z from state cs

global nq

R = 0.5*numel(z);
y = norm(sum([nq*evan(z(1:R),z(R+1:2*R)),-cs],2));

end

function [y,yeq] = constraint(z)
% positive when any component of z has norm less than 0.1

global nq

R = 0.5*numel(z);
ncpts = sum(abs(nq*evan(z(1:R),z(R+1:2*R))).^2);
y = 0.1^2 - min(ncpts);
yeq = [];

end

end