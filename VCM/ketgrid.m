% VCM5  orthonomal basis solver

global N;  N = 70;  brackets

% Peter's Hamiltonian
% nhn = -4*(0:N)' gives an interesting sudden divergence

nhn = nhn/2;
nhn = nhn-4*(0:N)';

T=2*pi;  h=0.005;		% time axis
T = 0.1;
t = h*(0:ceil(T/h));
snapshots = h*(0:5);

iters = 4;                                          %%iterations of ODE solver

q0 = nq*evan(a0,'even');	% expansion of |a0> over number states

% initialise storage for data to plot

c = nan(N+1, length(t));
cc = nan(N+1, length(t), iters);

BUF = nan(1,length(t));
alpha.o = BUF;  number = BUF;  csize = BUF;  qsize = BUF;
rsdl = BUF;  rsdln = BUF;

% initial condition

c(:,1) = q0;

% Exact values for comparison

qe = exp(-1i*nhn*t).*repmat(q0,size(t));	% column i is q(t(i))
alpha.e = sum(conj(qe).*(aop*qe));

% Integration loop

for i = 1:length(t) 

	% collect data to plot

	qo = c(:,i);
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
		dc = -1i*h*nhn.*ch;
	end
    	c(:,i+1) = c(:,i) + dc;

end

% Eigenvalues of discretised Hamiltonian and step operators

M = diag(1i*nhn);  ew = eig(M);
figure, plot(1:(N+1),real(ew),'+k',1:(N+1),imag(ew),'.k',1:(N+1),nhn,'or')
title 'Discretised H eigenvalues'

P = 1+h*M;  ew = flip(eig(P));
figure, semilogy(1:(N+1),sort(abs(ew)),'.k'), hold on
title 'size of semi-implicit step eigenvalues'

P = 1+h*M*P;  ew = flip(eig(P));
semilogy(1:(N+1),sort(abs(ew)),'.b')

P = 1+h*M*P;  ew = flip(eig(P));
semilogy(1:(N+1),sort(abs(ew)),'.g')

P = 1+h*M*P;  ew = flip(eig(P));
semilogy(1:(N+1),sort(abs(ew)),'.r')
legend('euler', '2', '3', '4 itns', 'Location', 'SouthEast')


% set up plotting grid

x = -10:0.3:10;  y = -10:0.3:10;
[X,Y] = meshgrid(x,y);  Z = X(:)+1i*Y(:);
Aps = nq*evan(Z,'even');

% plot derivatives

figure, title derivatives
for ti = snapshots, for j = 1:iters

	[~,i] = min(abs(t-ti));
	subplot(length(snapshots), iters, iters*(i-1) + j)
	zplot(x,y,Aps'*(cc(:,i+1,j)-c(:,i))), axis equal
	text(-4.4,8,num2str(norm(cc(:,i+1,j)-c(:,i))), 'Color', 'white')
	if j == 1
		ylabel(sprintf('t=%.3f, |c|=%.1e', t(i), norm(c(:,i))))
	end

end, end

 
% plot snapshots in phase space

% for ti = snapshots
% 
% 	[~,i] = min(abs(t-ti));
% 	qo = A*c(:,i);
% 
% 	figure, subplot 331
% 	plog(x,y,Aps'*qe(:,i),a), title(sprintf('exact state t = %.1e', t(i)))
% 
% 	subplot 332, plog(x,y,Aps'*A*pinvA*qe(:,i),a)
% 	text(-4.4,8,sprintf('expansion size %.1e', norm(pinvA*qe(:,i))), 'Color', 'white')
% 	title 'frame expansion'
% 	
% 	subplot 333, plog(x,y,Aps'*(qe(:,i)-A*pinvA*qe(:,i)),a)
% 	text(-4.4,8,sprintf('residual %.1e', norm(q0-A*c(:,1))), 'Color', 'white')
% 	title 'residual'
% 	
% 	subplot 334, plog(x,y,Aps'*(nhn.*qe(:,i)),a), title 'exact derivative'
% 	text(-4.4,8,sprintf('derivative size %.1e', norm(nhn.*qe(:,i))), 'Color', 'white')
% 	title 'exact derivative'
% 
% 	subplot 335, plog(x,y,Aps'*A*pinvA*(nhn.*qe(:,i)),a)
% 	text(-4.4,8,sprintf('expansion size %.1e', norm(pinvA*(nhn.*qe(:,i)))), 'Color', 'white')
% 	
% 	subplot 336, plog(x,y,Aps'*(nhn.*qe(:,i)-A*pinvA*(nhn.*qe(:,i))),a)
% 	text(-4.4,8,sprintf('residual %.1e', norm(nhn.*qe(:,i)-A*pinvA*(nhn.*qe(:,i)))), 'Color', 'white')
% 	
% 	subplot 337, plog(x,y,Aps'*qo,a), title 'simulated state'
% 	text(-4.4,8,sprintf('expansion size %.1e', norm(c(:,i))), 'Color', 'white')
% 
% 	subplot 338, plog(x,y,Aps'*(nhn.*qo),a)
% 	title 'simulated derivative'
% 	
% 	subplot 339, plog(x,y,Aps'*(qe(:,i)-qo),a)
% 	text(-4.4,8,sprintf('residual %.1e', norm(qe(:,i)-qo)), 'Color', 'white')
% 	title 'simulated state residual'
% 
% end

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

ccc = permute(cc, [1 3 2]);  ccc = reshape(ccc, N+1, []);
figure, waterfall(log(abs(ccc))')

%figure, fi = gcf;
%lines = semilogy(t, component_lengths);
%set(lines, 'Color', 'k');
%set(lines, {'LineStyle'}, {'-' '--' '--' '-'}');
%xlabel t
%ylabel lengths
%legend variational gabor