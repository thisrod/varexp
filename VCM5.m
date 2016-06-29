% VCM5  Coherent state Dirac-Frenkel solver for the quartic oscillator
%
% By Peter Drummond, ported to Matlab by Rodney Polkinghorne

% Independent parameters

% N.B: nhn = -4*(0:N)' gives an interesting sudden divergence

warning off MATLAB:rankDeficientMatrix

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
c = 0;		% get rid of this once loop is refactored

T=2*pi;  h=0.005;		% time axis
t = h*(0:ceil(T/h));

iters=4;                                          %%iterations of ODE solver

% initialise storage for data to plot

one.o = nan(1,length(t));  one.c =  one.o;  one.p = one.o;	% Observed and Comparision
two = one;  three = one;  four = one;
zdiscrep = one.o; adiscrep = one.o; hdiscrep = one.o; qdiscrep = one.o;  
csize = one;  qsize = one;
orank = one.o;  rrank = one.o;


zp = [zeros(R,1);  a0 + sigma*randn(R,2)*[1; 1i] ];	% draw initial ensemble

%% Integration loop

for i = 1:length(t) 

	% independent brackets
	qo = sum(nq*evan(zp), 2);
	qsize.o(i) = norm(qo);
	alphao = qo'*aop*qo/qsize.o(i)^2;
	one.o(i) = real(alphao);  two.o(i) = imag(alphao);
	three.o(i) = sum(abs(qo).^2.*(0:N)')/qsize.o(i)^2;
	csize.o(i) = norm(c);
  
	% Exact values for comparison
	
	qe = exp(-1i*nhn*t(i)).*c0;
	alphae = qe'*aop*qe;
	one.c(i)=real(alphae);
	two.c(i)=imag(alphae);
	three.c(i)=real(n0);
	
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
	orank(i) = rank(Dq);   rrank(i) = rank([Dq; sqrt(epsilon)*eye(2*R)]);  

end


figure
plot(t,orank,':k',t,rrank,'-k');
xlabel t
ylabel rank
legend original regularised

figure
plot(t,one.o,'-r', t,one.c,'-k', t, one.p, ':k');
xlabel t
ylabel X_1

figure
plot(t,two.o,'-r',t,two.c,'-k', t, two.p, ':k');            %%Plot quadrature Y
xlabel t
ylabel X_2

figure
plot(t,three.o,'-r', t,three.c,'-k', t,three.p,':k');             %%Plot mean N
xlabel t
ylabel N

figure
plot(t,qsize.o,'-r', t,qsize.p,':k');
xlabel t
ylabel 'norm of |\psi>'

figure
plot(t,csize.o,'-r', t,csize.p,':k');                           %%Plot max_W
xlabel t
ylabel '|c|'

