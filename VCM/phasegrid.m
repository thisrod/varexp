% VCM5  Coherent state grid solver for the quartic oscillator

global N;  N = 70;  brackets

% Peter's Hamiltonian
% nhn = -4*(0:N)' gives an interesting sudden divergence

nhn = nhn/2;
nhn = nhn-4*(0:N)';

a0 = 2;		% coherent amplitude of initial state
% q0 = zeros(N+1,1);  q0(5) = 1;		% number state 2
q0 = nq*evan(a0,'even');	% expansion of |a0> over number states
epsilon=1i*1.e-4;                                 %%stabilize matrix

T=2*pi;  h=0.005;		% time axis
T = 0.05;
t = h*(0:ceil(T/h));
snapshots = h*(0:5);

iters = 4;                                          %%iterations of ODE solver

% set up coherent state grid

l = 1.1;			% frame density, a0 is not a grid point
[ax,ay] = meshgrid(-4:l:4);  a = ax(:)+1i*ay(:);
a = a(abs(a) <= 4);
A = nq*evan(a,'even');
pinvA = pinv(A);

% initialise storage for data to plot

c = nan(length(a), length(t));
cc = nan(length(a), length(t), iters);

BUF = nan(1,length(t));
alpha.o = BUF;  number = BUF;  csize = BUF;  qsize = BUF;
rsdl = BUF;  rsdln = BUF;
%urank = BUF;  rrank = BUF;
%fcondition = BUF;  vcondition = BUF;

% initial condition, the least squares expansion of |a0> over A

c(:,1) = pinvA*q0;

% Exact values for comparison

qe = exp(-1i*nhn*t).*repmat(A*c(:,1),size(t));	% column i is q(t(i))
alpha.e = sum(conj(qe).*(aop*qe))/norm(A*c(:,1))^2;

% Integration loop

for i = 1:length(t) 

	% collect data to plot

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
		dc = pinvA*(-1i*h*nhn.*(A*ch));
	end
    	c(:,i+1) = c(:,i) + dc;

end

% Eigenvalues of discretised Hamiltonian and step operators

M = pinvA*diag(-1i*nhn)*A;
figure, plot(1:length(a), sort(abs(eig(M))), 'vk', 1:length(a), sort(abs(nhn(1:length(a)))), '^k')
hold on, plot([1 length(a)], 2/h*[1 1], '-k')
title 'H eigenvalues and Nyquist limit'
legend('discrete', 'exact', 'Location', 'SouthEast')


% set up plotting grid

x = -10:0.3:10;  y = -10:0.3:10;
[X,Y] = meshgrid(x,y);  Z = X(:)+1i*Y(:);
Aps = nq*evan(Z,'even');

% Maximum eigenvector of discrete H

[ev,ew] = eig(M);  ew = diag(ew);
[~,i] = max(abs(ew));
figure, plog(x,y,Aps'*A*ev(:,i),a)
title 'largest ev of discrete H'

figure, semilogy(0:N, abs(A*ev(:,i))/norm(A*ev(:,i)), '.k')
title 'number components of largest ev'


% plot derivatives

figure, title derivatives
for ti = snapshots, for j = 1:iters

	[~,i] = min(abs(t-ti));
	subplot(length(snapshots), iters, iters*(i-1) + j)
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

	figure, subplot 331
	plog(x,y,Aps'*qe(:,i),a), title(sprintf('exact state t = %.1e', t(i)))

	subplot 332, plog(x,y,Aps'*A*pinvA*qe(:,i),a)
	text(-4.4,8,sprintf('expansion size %.1e', norm(pinvA*qe(:,i))), 'Color', 'white')
	title 'frame expansion'
	
	subplot 333, plog(x,y,Aps'*(qe(:,i)-A*pinvA*qe(:,i)),a)
	text(-4.4,8,sprintf('residual %.1e', norm(q0-A*c(:,1))), 'Color', 'white')
	title 'residual'
	
	subplot 334, plog(x,y,Aps'*(nhn.*qe(:,i)),a), title 'exact derivative'
	text(-4.4,8,sprintf('derivative size %.1e', norm(nhn.*qe(:,i))), 'Color', 'white')
	title 'exact derivative'

	subplot 335, plog(x,y,Aps'*A*pinvA*(nhn.*qe(:,i)),a)
	text(-4.4,8,sprintf('expansion size %.1e', norm(pinvA*(nhn.*qe(:,i)))), 'Color', 'white')
	
	subplot 336, plog(x,y,Aps'*(nhn.*qe(:,i)-A*pinvA*(nhn.*qe(:,i))),a)
	text(-4.4,8,sprintf('residual %.1e', norm(nhn.*qe(:,i)-A*pinvA*(nhn.*qe(:,i)))), 'Color', 'white')
	
	subplot 337, plog(x,y,Aps'*qo,a), title 'simulated state'
	text(-4.4,8,sprintf('expansion size %.1e', norm(c(:,i))), 'Color', 'white')

	subplot 338, plog(x,y,Aps'*(nhn.*qo),a)
	title 'simulated derivative'
	
	subplot 339, plog(x,y,Aps'*(qe(:,i)-qo),a)
	text(-4.4,8,sprintf('residual %.1e', norm(qe(:,i)-qo)), 'Color', 'white')
	title 'simulated state residual'

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

%figure
%plot(t, 2*R-urank, ':k', t, 2*R-rrank, '-k');
%xlabel t
%ylabel 'rank deficiency in |\psi''>'
%legend original regularised

%figure
%semilogy(t, fcondition, '-k', t, vcondition, ':k');
%xlabel t
%ylabel 'frame condition'
%legend Gabor variational

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
plot(t, real(alpha.o), '-r', t, real(alpha.e),' -k');
xlabel t
ylabel X_1

figure
plot(t, imag(alpha.o), '-r', t, imag(alpha.e),' -k');
xlabel t
ylabel X_2

figure
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