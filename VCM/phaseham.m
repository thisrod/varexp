%PHASEHAM analyse the Kerr Hamiltionian discretised with wave packets

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

% The Mev turn out identical for large n, and numerically better
% So we'll use those instead.
[Gev,~] = eig(G);
Gkets = A*Gev;  Gkets = Gkets / diag(norms(Gkets));
Gns = sum(diag(0:N)*abs(Gkets).^2);
[Gns, P] = sort(Gns);  Gkets = Gkets(:,P);
Gew = diag(Mev'*G*Mev);

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
errorbar(0:R-1, Mns, sqrt(Mvarn), 'k'), hold on
plot([0 R-1], [0 R-1], ':k')
ylabel '<n>',  xlabel 'eigenkets of discrete number operator'

figure
errorbar(0:R-1, Ans, sqrt(Avarn), 'k'), hold on
plot([0 R-1], [0 R-1], ':k')
ylabel '<n>',  xlabel 'left SVs of expansion operator'


eyeR = eye(R);
x2 = 0:0.1:R+10;  y2 = -3:0.1:3;
[X,Y] = meshgrid(x2,y2);  Z2 = X+1i*Y;

% figure	% from smm program 24
% plot(Mew, '.k'), hold on
% sigmin = nan(length(y2), length(x2));
% wb = waitbar(0, 'please wait ...');
% for j = 1:length(x2), waitbar(j/length(x2))
% 	for i = 1:length(y2), sigmin(i,j) = min(svd(Z2(i,j)*eyeR-M)); end
% end, close(wb)
% [C1, ctrs1] = contour(x2, y2, sigmin, 10.^(-6:0.5:0)); colormap([0 0 0])
% axis equal
% title 'pseudospectrum of discrete N'
% 
% figure	% from smm program 24
% plot(0:R-1, zeros(1,R), '.k'), hold on
% sigmin = nan(length(y2), length(x2));
% wb = waitbar(0, 'please wait ...');
% for j = 1:length(x2), waitbar(j/length(x2))
% 	for i = 1:length(y2), sigmin(i,j) = min(svd(Z2(i,j)*eyeR-diag(0:R-1))); end
% end, close(wb)
% [C2, ctrs2] = contour(x2, y2, sigmin, 10.^(-6:0.5:0)); colormap([0 0 0])
% axis equal
% title 'pseudospectrum of N'

xx = -200:1:1000;

sigmin = nan(length(xx), 2);	% from smm program 24
for j = 1:length(xx)
		sigmin(j,1) = min(svd(xx(j)*eyeR-G));
		sigmin(j,2) = min(svd(xx(j)*eyeR-diag(nhn(1:R))));
end

% figure
% plot(Gew, '.k'), hold on
% [C3, ctrs3] = contour(x3, y3, sigmin, [10.^(-6:0.5:-3) 0.002*(1:10)]); colormap([0 0 0])
% axis equal
% title ''

figure
semilogy(xx, sigmin(:,1), '-k', xx, sigmin(:,2), ':k')
xlabel 'real \lambda',  legend discrete exact Location SouthEast
title 'pseudospectra of hamiltonian'


return

figure
zplot(0:R-1, 0:R-1, Mkets'*Gkets), axis image
title 'brackets of discrete N and H eigenkets'

figure
zplot(0:R-1, 0:R-1, abs(Au(1:R,:)).^2), axis image
xlabel m, ylabel '|<n|u_m>|'
title 'Fock components of expansion operator LSVs'

figure
zplot(0:R-1, 0:R-1, abs(Mev'*Mev)), axis image
title 'inner products of discrete N ev'

figure
B = bar(0:R-1, abs(Mev'*Av).^2, 'stacked');
set(B, 'EdgeColor', 'w', 'FaceColor', [0.5 0.5 0.5])
set(B(R-1:R), 'FaceColor', [0.2 0.2 0.2])
xlim([-0.5 R-0.5]), ylim([0 1])
title 'expansions of discrete N evs over RSVs'

for i = [1:5 (R-5):R]
	figure, zplot(x,y,Aps'*A*Mev(:,i)), hold on, axis image
	title(['Discrete eigenstate number ' num2str(i-1)])
	text(-9,6, {sprintf('<n> = %.1f', Mns(i)); sprintf('\\Deltan = %.1f', sqrt(Mvarn(i)))}, ...
		'Color', 'white')

	figure, zplot(x,y,Aps'*Au(:,i)), hold on, axis image
	title(['Left singular ket number ' num2str(i-1)])
	text(-9,6, {sprintf('<n> = %.1f', Ans(i)); sprintf('\\Deltan = %.1f', sqrt(Avarn(i))); ...
		sprintf('\\sigma = %.1e', Asw(i))}, 'Color', 'white')
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
