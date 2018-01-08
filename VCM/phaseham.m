%PHASEHAM analyse the wave packet discretisation of the number operator

global N;  N = 70;  brackets

% Peter's Hamiltonian
nhn = nhn/2;
nhn = nhn-4*(0:N)';

laf = 3e-4;		% Lagrange multiplier for Tyc conditioning.
h = 0.005;		% time step

% phase space grid and expansion operator

l = 1.1;			% frame density
[ax,ay] = meshgrid(-4:l:4);  a = ax(:)+1i*ay(:);
a = a(abs(a) <= 4);
R = length(a);
[~, i] = sort(abs(a));  a = a(i);
A = nq*evan(a,'even');
pinvA = pinv(A);

% Discretised number operator and hamiltonian

M = pinvA*diag(0:N)*A;
G = pinvA*diag(nhn)*A;

% Eigenvalues are stored as vectors
% Eigenvectors are stored as columns, in order of Fock number
% Suffix convention: ev, ew, sw, u = lsv, v = rsv
%     c = expansion coefficients
% The discrete hamiltonian has eigenkets in Hilbert space, 
% and eigenvectors in C^R.  These are also the eigenkets of the
% discrete hamiltonian.

[Au,Asw,Av] = svd(A,'econ');  Asw = diag(Asw);

Hew = nhn(1:R);   % Hev is the unit matrix
Ans = sum(diag(0:N)*abs(Au).^2);
Avarn = sum(diag((0:N).^2)*abs(Au).^2) - Ans.^2;

[Mev,Mew] = eig(M);  Mew = diag(Mew);
Mkets = A*Mev;  Mkets = Mkets / diag(norms(Mkets));
Mns = sum(diag(0:N)*abs(Mkets).^2);
[Mns, P] = sort(Mns);
Mew = Mew(P);  Mev = Mev(:,P);  Mkets = Mkets(:,P);
Mvarn = sum(diag((0:N).^2)*abs(Mkets).^2) - Mns.^2;

nr = zeros(N+1,R);  nr(1:R,:) = diag(0:R-1);
fockc = pinvA*nr;  fockc = fockc / diag(norms(fockc));

% Tychonov regularised H
Hy = pinv([A; laf*eye(R)])*[diag(nhn)*A; zeros(R)];
Tew = diag(Mev'*Hy*Mev); % use norms

% phase space grid for plotting

x = -10:0.3:10;  y = -10:0.3:10;
[X,Y] = meshgrid(x,y);  Z = X(:)+1i*Y(:);
Aps = nq*evan(Z,'even');

ixs = [1:4 round(linspace(5,length(a),4))];
cmap = phase(128);
ff = linspace(0, 2*pi, size(cmap,1));

% Graphs, as described in titles.

figure
errorbar(0:R-1, Ans, sqrt(Avarn), 'k'), hold on
plot([0 R-1], [0 R-1], ':k'),  xlim([0 45])
title 'ranges of Fock kets in expansion singular kets'
ylabel 'n',  xlabel 'left singular kets of expansion operator'

figure
zplot(0:R-1, 0:R-1, abs(Au(1:R,:)).^2), axis image
xlabel m, ylabel '|<n|u_m>|'
title 'Fock components of expansion operator LSVs'

xx = -10:0.1:50;

sigmin = nan(length(xx), 2);	% from smm program 24
for j = 1:length(xx)
		sigmin(j,1) = min(svd(xx(j)*eye(R)-M));
		sigmin(j,2) = min(svd(xx(j)*eye(R)-diag(0:R-1)));
end

figure
semilogy(xx, sigmin(:,1), '-k', xx, sigmin(:,2), '-k'), ylim([1e-6 50])
xlabel 'real \lambda',  ylabel \epsilon
title 'pseudospectra of number operator and discretisation'

figure
zplot(0:R-1, 0:R-1, abs(Mev'*Mev)), axis image
title 'inner products of discrete N ev'

figure
B = bar(0:R-1, abs(fockc'*Av).^2, 'stacked');
set(B, 'EdgeColor', 'w', 'FaceColor', [0.1 0.1 0.1])
set(B(R-1:R), 'FaceColor', [0.5 0.5 0.5])
xlim([-0.5 R-0.5]), ylim([0 1])
title 'components of RSVs in least squares Fock kets'

return

for i = [1:5 (R-5):R]
	figure, zplot(x,y,Aps'*A*Mev(:,i)), hold on, axis image
	title(['Discrete eigenstate number ' num2str(i-1)])
	text(-9,6, {sprintf('<n> = %.1f', Mns(i)); sprintf('\\Deltan = %.1f', sqrt(Mvarn(i)))}, ...
		'Color', 'white')
end

% Singular vectors of A

for i = ixs
	figure, zplot(x,y,Aps'*Au(:,i)), hold on, axis image
	title(sprintf('Left singular ket number %d', i))
	text(-9,6, {sprintf('<n> = %.1f', Ans(i)); sprintf('\\Deltan = %.1f', sqrt(Avarn(i))); ...
		sprintf('\\sigma = %.1e', Asw(i))}, 'Color', 'white')
	scale = 10/max(abs(Av(:,i)));
	for j = 1:length(a)
		plot(a(j), 'o', 'MarkerSize', scale*abs(Av(j,i)), ...
			'MarkerEdgeColor', interp1(ff, cmap, mod(angle(Av(j,i)), 2*pi)))
	end
end

figure, subplot 311
plot(0:R-1, Gew, 'vk', 0:R-1, Hew, '^k')
hold on, plot([0 R-1], 2/h*[1 1], '-k')
text(1, 2/h, 'Nyquist limit')
title 'H eigenvalues'
legend('discrete', 'exact', 'Location', 'SouthEast')

subplot 312
semilogy(0:R-1, Asw, '.k', [0 R-1], [laf laf], '-k');
text(1, laf, '\lambda')
xlabel n, ylabel '\sigma_n'
title 'Singular values of expansion'

subplot 313
pew = Hew.*Asw.^2./(Asw.^2+laf.^2);
plot(0:R-1, Tew, 'vk', 0:R-1, Hew, '^k', 0:R-1, pew, '.k')
hold on, plot([0 R-1], 2/h*[1 1], '-k')
title 'regularised H eigenvalues (should compute from Mev)'
legend('discrete', 'exact', 'predicted', 'Location', 'NorthWest')
