% Expand the derivative of a coherent state over two samples, as they
% approach each other.

global N; N = 20; brackets

[~, f] = evan([0.1 -0.1], 'even');
f(1) = f(1)+1i*pi;
c = zeros(N+1,2);
c(:,1) = nq*evan(0,0);
c(:,2) = sum(nq*evan([0 1i*pi] ,[0.1 -0.1]), 2);
for j = 1:2; c(:,j) = c(:,j)/norm(c(:,j)); end

nh = 13;  h = logspace(-1.2,log10(2.5), nh);  
psaxis = -3:0.1:3;
[x,y] = meshgrid(psaxis);  grid = (nq*evan(x(:)+1i*y(:), 'even'))';

figure
nms = zeros(2,nh);  rds = nms;
for j = 1:2
for i=1:nh
    a = 0.5*[-h(i) h(i)];  [M,f0] = evan(a, 'even');
    d = (nq*M)\c(:,j);
    ch = sum(nq*evan(f0+log(d), a), 2);
    nms(j,i) = norm(d);  rds(j,i) = norm(c(:,j)-ch);
    np = 4*(i-1)/(nh-1);
    if mod(np, 1)==0    % Matlab's test for an integer
        subplot(4, 5, 11+5*(j-1)+np)
        zplot(psaxis, psaxis, grid*ch);
        hold on; plot([h(i) -h(i)], [0 0], 'ow')
        set(gca, 'DataAspectRatio', [1 1 1])
    end
end
end

subplot 411
plot(h, nms(1,:), '-k', h, nms(2,:), '-r')
ylabel('expansion norms')
subplot 412
plot(h, rds(1,:), '-k', h, rds(2,:), '-r')
ylabel('expansion residuals')
xlabel('h')

print(gcf, 'cancel.pdf', '-dpdf')