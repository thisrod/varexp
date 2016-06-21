% VCM5  Coherent state Dirac-Frenkel solver for the quartic oscillator
%
% By Peter Drummond, ported to Matlab by Rodney Polkinghorne

% Independent parameters

a0=2;		% coherent amplitude of initial state
epsilon=1i*1.e-4;                                 %%stabilize matrix
%%om= [0;0;0];					                  %%diagonal frequency
om= [0;-4];
%%kap=[0,1,1];                                      %%anharmonic coupling
kap=[0,1]; 
R = 8;		% variational components
isup = 1;		% index into stores - remove

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
%%itsteep=0;                                        %%steepest descent iterations
nmeasure=4;                                       %%number of measurements
store=zeros(length(t),nsup,nmeasure);
exact=zeros(length(t),nmeasure);
NM=R*M;                                           %%multi-vector size
omexp=eye(M,M);
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
nstore=zeros(length(t),R);
wstore=zeros(length(t),R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 	Starts loop in time:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

for it=1:length(t);                                      %%loop until time T
  if it>1                                         %%If first time, initia0ize
    a=omexp*a;
    a1=a;                                         %%Store multivector
    vec1=zeros(NM,1);
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
       a=a1+da;
    end;                                          %%end iterations
    a=omexp*(a1+2.*da); 
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
    store(it,isup,1)=store(it,isup,1)+sum(real(rho(m,:).*a(2,:)));	
    store(it,isup,2)=store(it,isup,2)+sum(real(rho(m,:).*a(2,:)/1i));
    store(it,isup,3)=store(it,isup,3)+real(sum(rho(m,:).*a(2,:))*a(2,m)');           	   
  end;
  store(it,isup,1)=store(it,isup,1)/rhosum;
  store(it,isup,2)=store(it,isup,2)/rhosum;
  store(it,isup,3)=store(it,isup,3)/rhosum;
  store(it,isup,4)=rhosum;
%  t(it)=(it-1)*h;                                %%Store t-coordinate
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%% 	Ca0culates exact observables for initia0 coherent state case; assumes
%%  a single anharmonic oscillator with nonlinear coupling and detuning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
ax=a0(1)*exp(n0*(exp(-1i*kap(2)*t(it))-1)-1i*om(2)*t(it));
exact(it,1)=real(ax);
exact(it,2)=imag(ax);
exact(it,3)=real(n0); 

end;                                              %%end time loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%%        Output section - gives numerica0 output data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%%        Graphics section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

fs=20;
figure(1);
plot(t,store(:,1,1),'--',t,store(:,2,1),'-',t,exact(:,1),'-.');          %%Plot quadrature X
set(gca,'FontSize',fs);
xlabel('t ','FontSize',fs);
ylabel('X','FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

figure(2);
plot(t,store(:,1,2),'--',t,store(:,2,2),'-',t,exact(:,2),'-.');            %%Plot quadrature Y
set(gca,'FontSize',fs);
xlabel('t ','FontSize',fs);
ylabel('Y','FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

figure(3);
plot(t,store(:,1,3),'--',t,store(:,2,3),'-',t,exact(:,3),'-.');             %%Plot mean N
set(gca,'FontSize',fs);
xlabel('t ','FontSize',fs);
ylabel('N','FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

figure(4);
set(gca,'FontSize',fs);
plot(t,store(:,1,4),'--',t,store(:,2,4),'-');                           %%Plot norm
xlabel('t ','FontSize',fs);
ylabel('Rho','FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

figure(5);
set(gca,'FontSize',fs);
plot(t,nstore);                           %%Plot max_N
xlabel('t ','FontSize',fs);
ylabel('N_values','FontSize',fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

figure(6);
set(gca,'FontSize',fs);
plot(t,wstore);                           %%Plot max_W
xlabel('t ','FontSize',fs);
ylabel('Weights','FontSize',fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

