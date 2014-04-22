% Compare approximation by 2*R components to R components and
% their derivatives

global N;  N = 30;  brackets
R = 20;  a = 2*sqrt(0.5)*[1 1i]*randn(2,R);  gabor
mesh = -5:0.2:5;  phasespace

strides = [4 2 2 1];
expansions = {DA(:,1:4:end) A(:,1:2:end) DA(:,1:2:end) A};

figure
ens = a;
plot(real(ens), imag(ens), 'ok');  hold on
ens = a(1:2:end);
plot(real(ens), imag(ens), 'ok', 'MarkerFaceColor', [0.7 0.7 0.7])
ens = a(1:4:end);
plot(real(ens), imag(ens), 'ok', 'MarkerFaceColor', 'k')
set(gca, 'DataAspectRatio', [1 1 1])


figure
for i=1:4
    subplot(2,2,i)
    M = expansions{i};  zplot(0:N, 0:N, eye(N+1)-M*pinv(M))
    switch i
        case 1
            title('Expansion')
            ylabel([num2str(0.5*R) ' vectors'])
        case 2
            title('Variation')
        case 3
            ylabel([num2str(R) ' vectors'])
    end
end

figure
for i=1:4
    [U,S,V] = svd(expansions{i}, 'econ');
    subplot(4,R+1,R*i+i-R)
    ns = numel(diag(S));
    semilogy(0:ns-1, diag(S), '+k')
    for j=1:ns
        subplot(4,R+1,R*i+i-R+j)
        zplot(mesh, mesh, Aps'*U(:,j)); hold on
        plot(real(a(1:strides(i):end)), imag(a(1:strides(i):end)), '.w')
        set(gca,'DataAspectRatio',[1 1 1])
    end
end