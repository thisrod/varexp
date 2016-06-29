% VCM5  Coherent state Dirac-Frenkel solver for the quartic oscillator
%
% By Peter Drummond, ported to Matlab by Rodney Polkinghorne


global N;  N = 50;  brackets

% Peter's Hamiltonian
% nhn = -4*(0:N)' gives an interesting sudden divergence

nhn = nhn/2;
nhn = nhn-4*(0:N)';

a0=2;		% coherent amplitude of initial state
epsilon=1i*1.e-4;                                 %%stabilize matrix

R = 16;		% variational components

% Derived parameters

n0 = abs(a0)^2;

sigma=.02;                                        %%initial standard deviation

T=2*pi;  h=0.005;		% time axis
t = h*(0:ceil(T/h));
snapshots = [0 0.1 T/2];

iters=4;                                          %%iterations of ODE solver

% initialise storage for data to plot

z = nan(2*R, length(t));

BUF = nan(1,length(t));
alpha.o = BUF;  number = BUF;  csize = BUF;  qsize = BUF;
urank = BUF;  rrank = BUF;
fcondition = BUF;  vcondition = BUF;

% draw initial ensemble

z(:,1) = [zeros(R,1);  a0 + sigma*randn(R,2)*[1; 1i] ];
c0 = sum(nq*evan(z(:,1)), 2);

% Integration loop

for i = 1:length(t) 

	% independent brackets
	qo = sum(nq*evan(z(:,i)), 2);
	c = sqrt(sum(abs(nq*evan(z(:,i))).^2));
	qsize(i) = norm(qo);
	alpha.o(i) = qo'*aop*qo/qsize(i)^2;
	number(i) = sum(abs(qo).^2.*(0:N)')/qsize(i)^2;
	csize(i) = norm(c);
 	
	if i == length(t), break, end
	
	zh = z(:,i);
	dz = zeros(2*R,1);
	for j = 1:iters
		c = sqrt(sum(abs(nq*evan(zh)).^2));
		w = zh;  w(1:R) = w(1:R) - log(max(c));
		ensemble = nq*evan(w);
		Dq = [ensemble, ndqr*evan(w)];
		AA = Dq'*Dq;
		Hq = nhn.*sum(ensemble,2);
		HH = Dq'*Hq;
%		dz = dz + [Dq; sqrt(epsilon)*eye(2*R)] \ [-1i*Hq*h/2-Dq*dz; zeros(2*R,1)];
		dz = dz + (AA+epsilon*eye(2*R))\(-1i*HH*h-AA*dz);
		zh = z(:,i) + dz/2;
	end
    	z(:,i+1) = z(:,i) + dz;
	urank(i) = rank(Dq);   rrank(i) = rank([Dq; sqrt(epsilon)*eye(2*R)]);
	fcondition(i) = cond(ensemble);  vcondition(i) = cond(Dq);

end
 
% Exact values for comparison

qe = exp(-1i*nhn*t).*repmat(c0,size(t));	% column i is q(t(i))
alpha.e = sum(conj(qe).*(aop*qe))/norm(c0)^2;

% Set up phase space grid

x = -5:0.2:5;  y = -5:0.2:5;
[X,Y] = meshgrid(x,y);  Z = X(:)+1i*Y(:);
Aps = nq*evan(Z,'even');

% plot snapshots in phase space

for ti = snapshots

	[~,i] = min(abs(t-ti));
	ensemble = nq*evan(z(:,i));
	% find residuals by taking components in the orthogonal space
	[U,~,~] = svd(ensemble);  U = U(:,rank(ensemble)+1:end);
	rsdl = U'*Aps;  rsdl = sqrt(sum(abs(rsdl).^2));  rsdl = reshape(rsdl, size(X));
	nrms = pinv(ensemble)*Aps;  nrms = sqrt(sum(abs(nrms).^2));
	nrms = reshape(nrms, size(X));
	
	figure, zplot(x,y,Aps'*sum(ensemble, 2)), axis equal, hold on
	plot(z(R+1:end,i),'ow')
	contour(x,y,rsdl,[0.5 0.5],'-w')
	contour(x,y,nrms, csize(i)/qsize(i)*[1 1],'-g')
	contour(x,y,nrms, csize(i)/qsize(i)*[10 10],'-y')
	contour(x,y,nrms, csize(i)/qsize(i)*[100 100],'-r')
	title(sprintf('snapshot at t = %.2f', ti))

end

figure
plot(t, 2*R-urank, ':k', t, 2*R-rrank, '-k');
xlabel t
ylabel 'rank deficiency in |\psi''>'
legend original regularised

figure
semilogy(t, fcondition, '-k', t, vcondition, ':k');
xlabel t
ylabel 'frame condition'
legend Gabor variational

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
plot(t, csize, '-r');                           %%Plot max_W
xlabel t
ylabel '|c|'

