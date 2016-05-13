% VCM5  Coherent state Dirac-Frenkel solver for the quartic oscillator
%
% By Peter Drummond, ported to Matlab by Rodney Polkinghorne

% Independent parameters

a0=2;		% coherent amplitude of initial state
epsilon=1i*1.e-4;                                 %%stabilize matrix
%%om= [0;0;0];					                  %%diagona0 frequency
om= [0;-4];
%%kap=[0,1,1];                                      %%anharmonic coupling
kap=[0,1]; 

% Derived parameters

n0 = abs(a0)^2;
av=conj([0,a0]');                               %%combined initia0 vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%%	Sets up initia0 integration data va0ues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

dt=0.005;                                        %%time step 
sigma=.02;                                        %%initia0 standard deviation
tmax=2.*pi;                                       %%maximum time 
%%tmax=.1;
nt=floor(1+tmax/dt);                              %%points in time for dynamics
dt=tmax/(nt-1);                                   %%corrected time step 
modes=1;                                          %%number of modes 
M=1+modes;                                        %%size of state vector 
iters=4;                                          %%iterations of ODE solver
%%itsteep=0;                                        %%steepest descent iterations
nsup=2;                                           %%number of superposition loops
nmeasure=4;                                       %%number of measurements
store=zeros(nt,nsup,nmeasure);
t=zeros(nt,1);
exact=zeros(nt,nmeasure);
for isup=1:nsup;                                  %%loop over superpositions
N=8*isup;                                         %%number of superpositions
NM=N*M;                                           %%multi-vector size
omexp=eye(M,M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%%	Sets up integration vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

da=0.0*1i+zeros(M,N);                             %%complex matrix of zeros 
a = repmat(av,1,N) + sigma*(randn(M,N)+1i*randn(M,N));                          %%Initia0ises amplitudes
phimask=[ones(1),zeros(modes,1)']'*ones(1,N);
amask=[zeros(1),ones(modes,1)']'*ones(1,N);
a=a.*amask;                                       %%Mask phase to zero
a_0=a;                                            %%Stores initia0 amplitude
delt=eye(M,M);
deltN=eye(N,N);
delt(1,1)=0.0;
hd=ones(1,M);
hd(1)=0;
A1=zeros(NM,M,N)+1i;                                %%Initia0ises A matrix 
A=zeros(NM,NM)+1i;                                %%Initia0ises A matrix 
H=zeros(NM,1)+1i;                                 %%Initia0ises H matrix 
AD=zeros(NM,NM)+1i;                                %%Initia0ises A matrix 
HD=zeros(NM,1)+1i;                                 %%Initia0ises H matrix 
eps=epsilon*eye(NM,NM); 
nstore=zeros(nt,N);
wstore=zeros(nt,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 	Starts loop in time:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

for it=1:nt;                                      %%loop until time tmax
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
        rhod=rho.*deltN;
        for m=1:N;                                %%Count superpositions!
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
        vec1=vec1+(A+eps)\(-1i*H*dt/2-A*vec1); 
       da=reshape(vec1,M,N);                      %%Reshape to matrix     
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
  for m=1:N;                                      %%Count superpositions!
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
  t(it)=(it-1)*dt;                                %%Store t-coordinate
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%% 	Ca0culates exact observables for initia0 coherent state case; assumes
%%  a single anharmonic oscillator with nonlinear coupling and detuning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
ax=a0(1)*exp(n0*(exp(-1i*kap(2)*t(it))-1)-1i*om(2)*t(it));
exact(it,1)=real(ax);
exact(it,2)=imag(ax);
exact(it,3)=real(n0); 

end;                                              %%end time loop
end;						                      %%end superposition loop
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
ylabel('N_va0ues','FontSize',fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

figure(6);
set(gca,'FontSize',fs);
plot(t,wstore);                           %%Plot max_W
xlabel('t ','FontSize',fs);
ylabel('Weights','FontSize',fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

