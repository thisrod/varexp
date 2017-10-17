% VCM5  Coherent state grid solver for the quartic oscillator

global N;  N = 70;  brackets

% Peter's Hamiltonian
% nhn = -4*(0:N)' gives an interesting sudden divergence

nhn = nhn/2;
nhn = nhn-4*(0:N)';

a0 = 2;		% coherent amplitude of initial state
% q0 = zeros(N+1,1);  q0(5) = 1;		% number state 2
q0 = nq*evan(a0,'even');	% expansion of |a0> over number states

T=2*pi;  h=0.005;		% time axis
T = 0.05;
t = h*(0:ceil(T/h));
snapshots = [h 2*h 0.1 1 2 4];

gfac = 1.5;		% relative coefficient growth limit
% gfac = inf;

laf = 3e-4;		% Lagrange multiplier for Tyc conditioning.
% laf = 3e-4;	% value where parasitic solution takes over at t=pi
% laf = 0;

iters = 4;                                          %%iterations of ODE solver

% set up coherent state grid

l = 1.1;			% frame density, a0 is not a grid point
[ax,ay] = meshgrid(-4:l:4);  a = ax(:)+1i*ay(:);
a = a(abs(a) <= 4);
[~, i] = sort(abs(a));  a = a(i);
A = nq*evan(a,'even');
pinvA = pinv(A);

% initialise storage for data to plot
% cc(:,i,:) converge to c(:,i)

c = nan(length(a), length(t));
cc = nan(length(a), length(t), iters);

BUF = nan(1,length(t));
alpha.o = BUF;  number = BUF;  csize = BUF;  qsize = BUF;
rsdl = BUF;  rsdln = BUF;

% initial condition, the least squares expansion of |a0> over A

c(:,1) = pinvA*q0;
glim = gfac*norm(c(:,1));

% Exact values for comparison

qe = exp(-1i*nhn*t).*repmat(A*c(:,1),size(t));	% column i is q(t(i))
alpha.e = sum(conj(qe).*(aop*qe))/norm(A*c(:,1))^2;

% Integration loop

for i = 1:length(t) 

	% collect state data

	qo = A*c(:,i);
	rsdl(i) = norm(qe(:,i)-qo);
	rsdln(i) = sum(abs(qe(:,i)-qo).^2.*(0:N)')/rsdl(i)^2;
	qsize(i) = norm(qo);
	alpha.o(i) = qo'*aop*qo/qsize(i)^2;
	number(i) = sum(abs(qo).^2.*(0:N)')/qsize(i)^2;
	csize(i) = norm(c(:,i));
 	
	if i == length(t), break, end

	% propagate state
	
	dc = 0;
	for j = 1:iters
		cc(:,i+1,j) = c(:,i) + dc;
		ch = c(:,i) + dc/2;
		% Tychonov
		dc = h*pinv([A; laf*eye(length(a))])*[diag(-1i*nhn)*A*ch; zeros(size(a))];
		if norm(c(:,i) + dc) > glim
			% N.B. dc never comes out orthogonal to c
			% solve real l.s. problem with constant |c|
			cri = [real(c(:,i)); imag(c(:,i))];
			% Householder reflect c = grad(|c|²) onto e1, drop first column
			U = eye(2*length(a));
			v = U(:,1);  if sign(cri(1)) < 0, v = -v; end
			v = v + cri/norm(cri);
			U = U - 2*v*v'/(v'*v);  U = U(:, 2:end);
			% solve constrained least squares problem
			rhs = -1i*h*nhn.*(A*ch);  rhs = [real(rhs); imag(rhs)];
			dc = U*pinv([real(A) -imag(A); imag(A) real(A)]*U)*rhs;
			dc = dc(1:length(a)) + 1i*dc(length(a)+1:end);
			% keyboard
		end
		
	end
    	c(:,i+1) = c(:,i) + dc;

end

% Suffix convention: ev, ew, sw, u = lsv, v = rsv

Hew = nhn(1:length(a));
M = pinv(A)*diag(nhn)*A;
[ev,ew] = eig(M);  ew = diag(ew);
% ns = (norms(aop*A*ev)./norms(A*ev)).^2;
ns = sum(diag(0:N)*abs(A*ev).^2)./norms(A*ev).^2;
[ns, P] = sort(ns);
ev = ev(:,P);  ew = ew(P);
[Au,Asw,Av] = svd(A,'econ');  Asw = diag(Asw);
% Tychonov regularised H
Hy = pinv([A; laf*eye(length(a))])*[diag(nhn)*A; zeros(length(a))];
[Hyev, Hyew] = eig(Hy);  Hyew = diag(Hyew);
% yns = (norms(aop*A*Hyev)./norms(A*Hyev)).^2;
yns = sum(diag(0:N)*abs(A*Hyev).^2)./norms(A*Hyev).^2;
[yns, P] = sort(yns);
Hyev = Hyev(:,P);  Hyew = Hyew(P);

% Eigenvalues of discretised Hamiltonian and step operators

figure
zplot(1:length(a), 0:N, Au), axis image
title 'svs of A'

unit = eye(N+1,length(a));
Afock = pinv(A)*unit;
Afock = Afock/diag(norms(Afock));

figure, subplot 221
zplot(1:length(a), 1:length(a), ev), axis image
title 'evs of discretised H'

subplot 222
tmp2 = A*ev;  tmp2 = tmp2/diag(norms(tmp2));
imagesc(abs(tmp2(1:length(a),:)).^2), colormap gray, axis image
set(gca, 'YDir', 'normal')
title 'Fock discretised H'

subplot 223
zplot(1:length(a), 1:length(a), Afock), axis image
title 'number states'

subplot 224
tmp1 = A*Afock;  tmp1 = tmp1/diag(norms(tmp1));
imagesc(abs(tmp1(1:length(a),:)).^2), colormap gray, axis image
set(gca, 'YDir', 'normal')
title 'Fock number states'


figure
plot(0:length(a)-1, yns, '.r', 0:length(a)-1, ns, '.k', [0 length(a)], [0 length(a)], ':k')
title 'Average number in numerical eigenstates'
ylabel '<n>' , title 'excitation ns of discrete evs'
legend discrete regularised Location NorthWest

figure, subplot 311
plot(1:length(a), ew, 'vk', 1:length(a), Hew, '^k')
hold on, plot([1 length(a)], 2/h*[1 1], '-k')
text(1, 2/h, 'Nyquist limit')
title 'H eigenvalues'
legend('discrete', 'exact', 'Location', 'SouthEast')

subplot 312
semilogy(1:length(a), Asw, '.k', [0 length(a)], [laf laf]);
title 'Singular values of expansion'
xlabel n, ylabel '\sigma_n'

subplot 313
pew = Hew.*Asw.^2./(Asw.^2+laf.^2);

plot(1:length(a), Hyew, 'vk', 1:length(a), Hew, '^k', 1:length(a), pew, '.k')
hold on, plot([1 length(a)], 2/h*[1 1], '-k')
title 'regularised H eigenvalues'
legend('discrete', 'exact', 'predicted', 'Location', 'NorthWest')

% set up plotting grid

x = -10:0.3:10;  y = -10:0.3:10;
[X,Y] = meshgrid(x,y);  Z = X(:)+1i*Y(:);
Aps = nq*evan(Z,'even');

ixs = [1:4 round(linspace(5,length(a),4))];
cmap = phase(128);
ff = linspace(0, 2*pi, size(cmap,1));

% Singular vectors of A

for i = ixs
	figure, zplot(x,y,Aps'*Au(:,i)), hold on, axis image
	title(sprintf('Singular vector %d', i))
	text(-9,8,sprintf('\\sigma = %.1e', Asw(i)), 'Color', 'white')
	scale = 10/max(abs(Av(:,i)));
	for j = 1:length(a)
		plot(a(j), 'o', 'MarkerSize', scale*abs(Av(j,i)), ...
			'MarkerEdgeColor', interp1(ff, cmap, mod(angle(Av(j,i)), 2*pi)))
	end
end

% Eigenvectors of M

[~, P] = sort(abs(ew));
for i = ixs
	figure, zplot(x,y,Aps'*A*ev(:,P(i))), hold on, axis image
	title(sprintf('Discrete H eigenvector %d', i))
	text(-9,8,sprintf('ew = %.1e', ew(P(i))), 'Color', 'white')
	text(-9,7,sprintf('|A*ev| = %.1e', norm(A*ev(:,P(i)))), 'Color', 'white')
	scale = 10/max(abs(ev(:,P(i))));
	for j = 1:length(a)
		plot(a(j), 'o', 'MarkerSize', scale*abs(ev(j,P(i))), ...
			'MarkerEdgeColor', interp1(ff, cmap, mod(angle(ev(j,P(i))), 2*pi)))
	end
end

% expansions of number states

for i = ixs
	cs = pinv(A)*unit(:,i);
	figure, zplot(x,y,Aps'*A*cs), hold on, axis image
	title(sprintf('Number state %d', i))
	text(-9,7,sprintf('|c| = %.1e', norm(cs)), 'Color', 'white')
	scale = 10/max(abs(cs));
	for j = 1:length(a)
		plot(a(j), 'o', 'MarkerSize', scale*abs(cs(j)), ...
			'MarkerEdgeColor', interp1(ff, cmap, mod(angle(cs(j)), 2*pi)))
	end
end

% plot derivatives

figure, title derivatives
for ti = 1:length(snapshots), for j = 1:iters

	[~,i] = min(abs(t-snapshots(ti)));
	subplot(length(snapshots), iters, iters*(ti-1) + j)
	plog(x,y,Aps'*A*(cc(:,i+1,j)-c(:,i)),a)
	text(-4.4,8,num2str(norm(cc(:,i+1,j)-c(:,i))), 'Color', 'white')
	if j == 1
		ylabel(sprintf('t=%.3f, |c|=%.1e', t(i), norm(c(:,i))))
	end

end, end

 
% plot snapshots in phase space

for ti = snapshots

	[~,i] = min(abs(t-ti));
	qo = A*c(:,i);

	figure, subplot 241
	plog(x,y,Aps'*qe(:,i),a)
	text(-4.4,9,sprintf('state size %.1e', norm(qe(:,i))), 'Color', 'white')
	text(-4.4,7,sprintf('expansion size %.1e', norm(pinvA*qe(:,i))), 'Color', 'white')
	title(sprintf('state t = %.1e', t(i)))

	subplot 245, plog(x,y,Aps'*(qe(:,i)-A*pinvA*qe(:,i)),a)
	title(sprintf('residual fraction %.1e', norm(q0-A*c(:,1))/norm(qe(:,i))))
	
	subplot 242, plog(x,y,Aps'*(nhn.*qe(:,i)),a)
	text(-4.4,9,sprintf('derivative size %.1e', norm(nhn.*qe(:,i))), 'Color', 'white')
	text(-4.4,7,sprintf('expansion size %.1e', norm(pinvA*(nhn.*qe(:,i)))), 'Color', 'white')
	title 'derivative'
	
	subplot 246, plog(x,y,Aps'*(nhn.*qe(:,i)-A*pinvA*(nhn.*qe(:,i))),a)
	title(sprintf('%.1e', norm(nhn.*qe(:,i)-A*pinvA*(nhn.*qe(:,i)))/norm(nhn.*qe(:,i))))
	
	subplot 243, plog(x,y,Aps'*qo,a)
	text(-4.4,9,sprintf('state size %.1e', norm(A*c(:,i))), 'Color', 'white')
	text(-4.4,7,sprintf('expansion size %.1e', norm(c(:,i))), 'Color', 'white')
	title 'simulated state'
	
	subplot 247, plog(x,y,Aps'*(qe(:,i)-qo),a)
	title(sprintf('%.1e', norm(qe(:,i)-qo)/norm(qe(:,i))))

	subplot 244, plog(x,y,Aps'*(nhn.*qo),a)
	text(-4.4,9,sprintf('derivative size %.1e', norm(nhn.*qo)), 'Color', 'white')
	text(-4.4,7,sprintf('expansion size %.1e', norm(pinv(A)*(nhn.*qo))), 'Color', 'white')
	title 'simulated derivative'
	
	subplot 248, plog(x,y,Aps'*(nhn.*(qe(:,i)-qo)),a)
	title(sprintf('%.1e', norm(nhn.*(qe(:,i)-qo))/norm(nhn.*qe(:,i))))

%	% find residuals by taking components in the orthogonal space
%	[U,~,~] = svd(ensemble);  U = U(:,rank(ensemble)+1:end);
%	rsdl = U'*Aps;  rsdl = sqrt(sum(abs(rsdl).^2));  rsdl = reshape(rsdl, size(X));
%	nrms = pinv(ensemble)*Aps;  nrms = sqrt(sum(abs(nrms).^2));
%	nrms = reshape(nrms, size(X));
%	
%	figure, zplot(x,y,Aps'*sum(ensemble, 2)), axis equal, hold on
%	plot(z(R+1:end,i),'ow')
%	contour(x,y,rsdl,[0.5 0.5],'-w')
%	nrmax = component_lengths(3,i);
%	contour(x,y,nrms, nrmax*[1 1],'-g')
%	contour(x,y,nrms, nrmax*[10 10],'-y')
%	contour(x,y,nrms, nrmax*[100 100],'-r')
%	title(sprintf('state at t = %.2f', ti))
%	
%	D = [ensemble, ndqr*evan(z(:,i))];
%	[U,~,~] = svd(D);  U = U(:,rank(D)+1:end);
%	rsdl = U'*Aps;  rsdl = sqrt(sum(abs(rsdl).^2));  rsdl = reshape(rsdl, size(X));
%	nrms = pinv(D)*Aps;  nrms = sqrt(sum(abs(nrms).^2));
%	nrms = reshape(nrms, size(X));
%	
%	figure, zplot(x,y,Aps'*(nhn.*sum(ensemble, 2))), axis equal, hold on
%	plot(z(R+1:end,i),'ow')
%	contour(x,y,rsdl,[0.5 0.5],'-w')
%	nrmax = max(nrmax, component_lengths(4,i));
%	contour(x,y,nrms, nrmax*[1 1],'-g')
%	contour(x,y,nrms, nrmax*[10 10],'-y')
%	contour(x,y,nrms, nrmax*[100 100],'-r')
%	title(sprintf('derivative at t = %.2f', ti))

end

% plot vector sum diagrams (do this for snapshots)
    	
figure
for j = 1:length(snapshots)
	[~,i] = min(abs(t-snapshots(j)));
 	subplot(1, length(snapshots), j)
 	u = [c(:,i) c(:,i+1) squeeze(cc(:,i+1,:))];
 	[Q,R] = qr([real(u); imag(u)], 0);
 	plot([0 R(1,1)], [0 0], '-k', [R(1,1) R(1,2)], [0 R(2,2)], '-r', ...
 		R(1,3:end), R(2,3:end), 'xk', 'LineWidth', 2)
 	axis image
 	title(sprintf('t = %.4f', t(i)+h/2))
end

figure, subplot 311
semilogy(t, rsdl, '-k');
xlabel t
ylabel 'error ket size'

subplot 312
plot(t, rsdln, '-k');
xlabel t
ylabel '<n> for error ket'

subplot 313
semilogy(t, csize, '-k'), hold on
semilogy(t,norms(cc),'.k')
xlabel t
ylabel '|c|'

figure
subplot 311
plot(t, real(alpha.o), '-r', t, real(alpha.e),' -k');
xlabel t
ylabel X_1

subplot 312
plot(t, imag(alpha.o), '-r', t, imag(alpha.e),' -k');
xlabel t
ylabel X_2

subplot 313
plot(t, number, '-r', t([1 end]), number([1 1]), '-k');
xlabel t
ylabel <n>

figure
plot(t, qsize, '-r');
xlabel t
ylabel 'norm of |\psi>'


figure
plot(t, csize./qsize, '-r');
xlabel t
title 'ratio of coefficient norm to state norm'


ccc = permute(cc, [1 3 2]);  ccc = reshape(ccc, length(a), []);
figure, waterfall(0:N,h/4*(0:(size(ccc,2)-1)),log(abs(A*ccc))')
xlabel n, ylabel t
title 'Growth of parasitic eigenvector'

%figure, fi = gcf;
%lines = semilogy(t, component_lengths);
%set(lines, 'Color', 'k');
%set(lines, {'LineStyle'}, {'-' '--' '--' '-'}');
%xlabel t
%ylabel lengths
%legend variational gabor