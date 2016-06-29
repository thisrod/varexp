% VCM5  Coherent state Dirac-Frenkel solver for the quartic oscillator
%
% By Peter Drummond, ported to Matlab by Rodney Polkinghorne

% Independent parameters

% N.B: nhn = -4*(0:N)' gives an interesting sudden divergence

global N;  N = 30;  brackets
% Peter's Hamiltonian
nhn = nhn/2;
nhn = nhn-4*(0:N)';

a0=2;		% coherent amplitude of initial state
c0 = nq*evan(a0, 'even');	% expansion over number states
epsilon=1i*1.e-4;                                 %%stabilize matrix

R = 16;		% variational components

% Derived parameters

n0 = abs(a0)^2;

sigma=.02;                                        %%initial standard deviation
z = zeros(2*R, 1);

T=2*pi;  h=0.005;		% time axis
t = h*(0:ceil(T/h));

iters=4;                                          %%iterations of ODE solver

% initialise storage for data to plot

BUF = nan(1,length(t));
alpha.o = BUF;  number = BUF;  csize = BUF;  qsize = BUF;
urank = BUF;  rrank = BUF;

% draw initial ensemble
zp = [zeros(R,1);  a0 + sigma*randn(R,2)*[1; 1i] ];

% Integration loop

for i = 1:length(t) 

	% independent brackets
	qo = sum(nq*evan(zp), 2);
	c = sqrt(sum(abs(nq*evan(zp)).^2));
	qsize(i) = norm(qo);
	alpha.o(i) = qo'*aop*qo/qsize(i)^2;
	number(i) = sum(abs(qo).^2.*(0:N)')/qsize(i)^2;
	csize(i) = norm(c);
 	
	if i == length(t), break, end
	
	zph = zp;
	dzt = zeros(2*R,1);  dzp = dzt;	% Tychonov and Peter
	for j = 1:iters
		c = sqrt(sum(abs(nq*evan(zph)).^2));
		w = zph;  w(1:R) = w(1:R) - log(max(c));
		ensemble = nq*evan(w);
		Dq = [ensemble, ndqr*evan(w)];
		AA = Dq'*Dq;
		Hq = nhn.*sum(ensemble,2);
		HH = Dq'*Hq;
		dzt = dzt + [Dq; sqrt(epsilon)*eye(2*R)] \ [-1i*Hq*h/2-Dq*dzt; zeros(2*R,1)];
		dzp = dzp + (AA+epsilon*eye(2*R))\(-1i*HH*h-AA*dzp);
		zph = zp + dzp/2;
	end
    	zp = zp + dzp;
	urank(i) = rank(Dq);   rrank(i) = rank([Dq; sqrt(epsilon)*eye(2*R)]);  

end
 
% Exact values for comparison

qe = exp(-1i*nhn*t).*repmat(c0,size(t));	% column i is q(t(i))
alpha.e = sum(conj(qe).*(aop*qe));

figure
plot(t, 2*R-urank, ':k', t, 2*R-rrank, '-k');
xlabel t
ylabel 'rank deficiency in |\psi''>'
legend original regularised

figure
plot(t, real(alpha.o), '-r', t, real(alpha.e),' -k');
xlabel t
ylabel X_1

figure
plot(t, imag(alpha.o), '-r', t, imag(alpha.e),' -k');
xlabel t
ylabel X_2

figure
plot(t, number, '-r', t([1 end]), abs(a0)^2*[1 1]);
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

