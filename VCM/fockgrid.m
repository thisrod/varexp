% demonstrate parasitic solutions to the SI method

global N;  N = 25;  brackets

% Peter's Hamiltonian
nhn = nhn/2;
nhn = nhn-4*(0:N)';

a0 = 2;		% coherent amplitude of initial state
% q0 = zeros(N+1,1);  q0(5) = 1;		% number state 2
q0 = nq*evan(a0,'even');	% expansion of |a0> over number states

T=2*pi;  h=0.03;		% time axis
T = 0.09;
t = h*(0:ceil(T/h));
snapshots = [h 2*h 0.1 1 2 4];

iters = 4;                                          %%iterations of ODE solver

% initialise storage for data to plot
% cc(:,i,:) converge to c(:,i)

c = nan(N+1, length(t));
cc = nan(N+1, length(t), iters);

BUF = nan(1,length(t));
alpha.o = BUF;  number = BUF;  csize = BUF;  qsize = BUF;
rsdl = BUF;  rsdln = BUF;

% initial condition, the least squares expansion of |a0> over A

c(:,1) = q0;
cc(:,1) = q0;

% Exact values for comparison

qe = exp(-1i*nhn*t).*repmat(q0,size(t));	% column i is q(t(i))
alpha.e = sum(conj(qe).*(aop*qe))/norm(q0)^2;

% Integration loop

for i = 1:length(t) 

	% collect state data

	qo = c(:,i);
	rsdl(i) = norm(qe(:,i)-c(:,i));
	rsdln(i) = sum(abs(qe(:,i)-c(:,i)).^2.*(0:N)')/rsdl(i)^2;
	qsize(i) = norm(c(:,i));
	alpha.o(i) = qo'*aop*c(:,i)/qsize(i)^2;
	number(i) = sum(abs(c(:,i)).^2.*(0:N)')/qsize(i)^2;
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

% set(groot, 'defaultAxesFontName', 'Latin Modern Roman');
% set(groot, 'defaultAxesFontSize', 12);

% iterated hamiltonian ews

ww = linspace(0,3.5,50);
wm = 1-1i*ww;
for j = 2:4
	wm(j,:) = wm(j-1,:) + 2*(-1i*ww/2).^j;
end

figure
L = plot(ww, abs(wm), '-k');
[L.LineWidth] = deal(1,1.4,1.7,2);
xlabel '\it h\omega_k', ylabel '\it u_k^{\rm(\it j\rm)}'
ylim([0 6])
saveTightFigure ../wev.eps


% plot blowup

ccc = permute(cc, [1 3 2]);  ccc = reshape(ccc, N+1, []);
ccc = ccc(:,iters:end);	% trim nans that converge to initial value
ccc(:,1) = q0;

figure, S = surf(h/4*(0:(size(ccc,2)-1)),0:N,abs(ccc));
A = gca;  hold on
plot3(repmat([0 t(end)], N+1, 1), 0:N, abs(q0), '-k', 'LineWidth', 1.5);
psize = abs(q0).*(2*(h*nhn/2).^iters).^(t(end)/h);
plot3(repmat(t(end), N-16, 1), 17:N, ...
	psize(18:N+1), '--k', 'LineWidth', 1.5);
plot3([0 T], [17 17], q0(18)*[1 1], ':k', 'LineWidth', 1.5)
A.ZScale = 'log';  xlim([0 T]), view(-71.5,48)
 zlim([1e-6 10]),  A.XTick = 0:h:T,  A.ZTick = 10.^(-6:3:0)
% L = [A.YTickLabel; '\omega_n = 2/h'];  [A.YTick, i] = sort([A.YTick 17]);
% A.YTickLabel = L(i);
S.MeshStyle = 'row';  S.EdgeColor = 'k';  S.FaceColor = 'none';
xlabel '\it t', ylabel '\it n', zlabel '|\its_n\rm|'
% title 'Growth of parasitic eigenvectors'
saveTightFigure ../pev.eps

return

% set up plotting grid

x = -10:0.3:10;  y = -10:0.3:10;
[X,Y] = meshgrid(x,y);  Z = X(:)+1i*Y(:);
Aps = nq*evan(Z,'even');

% ixs = [1:4 round(linspace(5,R,4))];
cmap = phase(128);
ff = linspace(0, 2*pi, size(cmap,1));

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

