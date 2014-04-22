% Compute the frame bounds, and the dual frame, for a random ensemble of
% coherent states.

global N; N = 30; brackets
R = 15; a = 2*sqrt(0.5)*[1 1i]*randn(2,R); gabor
mesh = linspace(-5,5,20);  phasespace

[~, sv, ~] = svd(A);
S = A*A';
F = S\A;

figure
for r=1:R
    subplot(1,R+1,r+1)
    zplot(mesh, mesh, Aps'*F(:,r))
    set(gca, 'DataAspectRatio', [1 1 1])
end
subplot(1,R+1,1)
text(0.5,0.5, ['Frame bounds'; 'A = ' num2str(sv(end,end)); ...
    'B = ' num2str(sv(1,1))])