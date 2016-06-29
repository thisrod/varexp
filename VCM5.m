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
%%om= [0;0;0];					                  %%diagonal frequency
om= [0;-4];
%%kap=[0,1,1];                                      %%anharmonic coupling
kap=[0,1]; 
% kap = kap*0;
R = 16;		% variational components

% Derived parameters

n0 = abs(a0)^2;
av=conj([0,a0]');                               %%combined initial vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%%	Sets up initial integration data values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

sigma=.02;                                        %%initial standard deviation

T=2*pi;  h=0.005;		% time axis
t = h*(0:ceil(T/h));
modes=1;                                          %%number of modes 
M=1+modes;                                        %%size of state vector 
iters=4;                                          %%iterations of ODE solver

% data to plot
one.o = zeros(1,length(t));  one.c =  one.o;  one.p = one.o;	% Observed and Comparision
two = one;  three = one;  four = one;
zdiscrep = one.o; adiscrep = one.o; hdiscrep = one.o;

nstore =zeros(length(t),R);
wstore=zeros(length(t),R);

NM=R*M;                                           %%multi-vector size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%%	Sets up integration vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

da=0.0*1i+zeros(M,R);                             %%complex matrix of zeros 
a = repmat(av,1,R) + sigma*(randn(M,R)+1i*randn(M,R));                          %%Initia0ises amplitudes
phimask=[ones(1),zeros(modes,1)']'*ones(1,R);
amask=[zeros(1),ones(modes,1)']'*ones(1,R);
a=a.*amask;                                       %%Mask phase to zero
a_0=a;                                            %%Stores initial amplitude
delt=eye(M,M);
deltR=eye(R,R);
delt(1,1)=0.0;
hd=ones(1,M);
hd(1)=0;
A1=zeros(NM,M,R)+1i;                                %%Initia0ises A matrix 
A=zeros(NM,NM)+1i;                                %%Initia0ises A matrix 
H=zeros(NM,1)+1i;                                 %%Initia0ises H matrix 
AD=zeros(NM,NM)+1i;                                %%Initia0ises A matrix 
HD=zeros(NM,1)+1i;                                 %%Initia0ises H matrix 
eps=epsilon*eye(NM,NM); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 	Starts loop in time:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

for it=1:length(t);                                      %%loop until time T
  if it>1                                         %%If first time, initia0ize
    a1=a;                                         %%Store multivector
    vec1=zeros(NM,1);
    dzt = zeros(2*R,1);  dzp = dzt;	% Tychonov and Peter
    for iter=1:iters;                             %%iteration loop for ODE
        at=a.*amask+phimask;                      %%Store eigenva0ues
        at1=at+phimask;                           %%Store inner product terms
        logrho=0.5*(a'*at1+at1'*a);               %%Store log density matrix
        rhomax=max(max(real(logrho)));            %%Get maximum weight
        rho=exp(logrho-rhomax);                   %%Renorma0ize weight
        rhod=rho.*deltR;
        for m=1:R;                                %%Count superpositions!
          for k=1:M;                              %%Count modes!
            mk=(m-1)*M+k;
              for l=1:M;                          %%Count modes!	
                A1(mk,l,:)=rho(m,:).*(delt(k,l)+at(k,:)*at(l,m)');
              end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 	Ca0culate Hamiltonian for anharmonic oscillator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
   
h1=hd(k)*rho(m,:).*a(k,:).*(om(k)+kap(k)*a(k,m)'*a(k,:));            %%Commutator term
h2w= (a(:,m))'*diag(om)*a;	
h2k= (a(:,m).^2)'*diag(kap)*(a.^2);	
h2= rho(m,:).*at(k,:).*(h2w+0.5*h2k);                                        %%Hamiltonian term

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 	End problem-specific code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/ 	
              H(mk)=sum(h1+h2);
                  
          end;
        end;
        A=reshape(A1,NM,NM);
        %%B=A'; C=B*A; D=diag(B);
        vec1=vec1+(A+eps)\(-1i*H*h/2-A*vec1); 
       da=reshape(vec1,M,R);                      %%Reshape to matrix
       dza = da.';  dza = dza(:);
       
	z = a.';  z = z(:);
	c = sqrt(sum(abs(nq*evan(z)).^2));
	w = z;  w(1:R) = w(1:R) - log(max(c));
	ensemble = nq*evan(w);
	Dq = [ensemble, ndqr*evan(w)];
	AA = Dq'*Dq;
	Hq = nhn.*sum(ensemble,2);
	HH = Dq'*Hq;
	dzt = dzt + [Dq; sqrt(epsilon)*eye(2*R)] \ [-1i*Hq*h/2-Dq*dzt; zeros(2*R,1)];
	dzp = dzp + (AA+epsilon*eye(2*R))\(-1i*HH*h/2-AA*dzp);
       a=a1+da;
    end;                                          %%end iterations
    a=a1+2.*da;
    shuffle = [1:2:2*R 2:2:2*R];
    zdiscrep(it) = discrepency(dza, dzp);
    adiscrep(it) = discrepency(A(shuffle, shuffle), AA, 'fro');
    hdiscrep(it) = discrepency(H(shuffle), HH, 'fro');
    orank(it) = rank(Dq);   rrank(it) = rank([Dq; sqrt(epsilon)*eye(2*R)]);  
  end;                                            %%end if first time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%% 	Ca0culates observables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

  at1=a.*amask+2*phimask;                         %%Store eigenva0ues
  rho=0.5*(a'*at1+at1'*a);
  rho=exp(rho);                                   %%Renorma0ize weight
  rhosum=real(sum(sum(rho)));
  for m=1:R;                                      %%Count superpositions!
    nstore(it,m)=abs(a(2,m))^2;
    wstore(it,m)=real(rho(m,m));                                    
    one.o(it)=one.o(it)+sum(real(rho(m,:).*a(2,:)));	
    two.o(it)=two.o(it)+sum(real(rho(m,:).*a(2,:)/1i));
    three.o(it)=three.o(it)+real(sum(rho(m,:).*a(2,:))*a(2,m)');           	   
  end;
  one.o(it)=one.o(it)/rhosum;
  two.o(it)=two.o(it)/rhosum;
  three.o(it)=three.o(it)/rhosum;
  four.o(it)=rhosum;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%% 	Ca0culates exact observables for initia0 coherent state case; assumes
%%  a single anharmonic oscillator with nonlinear coupling and detuning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
ax=a0(1)*exp(n0*(exp(-1i*kap(2)*t(it))-1)-1i*om(2)*t(it));
ct = exp(-1i*nhn*t(it)).*c0;
alpha = ct'*aop*ct;
one.p(it)=real(ax);
one.c(it)=real(alpha);
two.p(it)=imag(ax);
two.c(it)=imag(alpha);
three.c(it)=real(n0); 

end;                                              %%end time loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%%        Output section - gives numerica0 output data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/


figure
plot(t,zdiscrep,'-k', t,adiscrep,':k',  t,hdiscrep,'--k');          %%Plot quadrature X
xlabel t, ylabel discrepency, legend z A H

figure
plot(t,orank,':k',t,rrank,'-k');          %%Plot quadrature X
xlabel t
ylabel rank
legend original regularised

figure
plot(t,one.o,'-r', t,one.c,'-k', t, one.p, ':k');          %%Plot quadrature X
xlabel t
ylabel X

figure
plot(t,two.o,'-r',t,two.c,'-k', t, two.p, ':k');            %%Plot quadrature Y
xlabel t
ylabel Y

figure
plot(t,three.o,'-r',t,three.c,'-k');             %%Plot mean N
xlabel t
ylabel N

figure
plot(t,four.o,'-k',t,four.c,'-r');                           %%Plot norm
xlabel t
ylabel Rho

figure
plot(t,nstore);                           %%Plot max_N
xlabel t
ylabel N

figure
plot(t,wstore);                           %%Plot max_W
xlabel t
ylabel Weights

