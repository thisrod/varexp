% sizes and residuals for expansions of coherent states

global N;  N = 60;  brackets
mesh = -5:0.2:5;  phasespace

R = 3;  h = 0.7;  [X,Y] = meshgrid(h*((0:R-1)-(R-1)/2));
a{1} = X(:) + 1i*Y(:);  A{1} = nq*evan(a{1}, 'even');
Z = R*h*(rand(R^2,2)-1/2);
a{2} = Z*[1; 1i];  A{2} = nq*evan(a{2}, 'even');
h = 1.5;  [X,Y] = meshgrid(h*((0:R-1)-(R-1)/2));
a{3} = X(:) + 1i*Y(:);  A{3} = nq*evan(a{3}, 'even');
Z = R*h*(rand(R^2,2)-1/2);
a{4} = Z*[1; 1i];  A{4} = nq*evan(a{4}, 'even');
R = 2;  h = 0.7*sqrt(2);  [X,Y] = meshgrid(h*((0:R-1)-(R-1)/2));
a{5} = X(:) + 1i*Y(:);  V = evan(a{5}, 'even');  A{5} = [nq*V, ndqr*V];
Z = R*h*(rand(R^2,2)-1/2);
a{6} = Z*[1; 1i];  V = evan(a{6}, 'even');  A{6} = [nq*V, ndqr*V];
R = 2;  h = 1.5*sqrt(2);  [X,Y] = meshgrid(h*((0:R-1)-(R-1)/2));
a{7} = X(:) + 1i*Y(:);  V = evan(a{7}, 'even');  A{7} = [nq*V, ndqr*V];
Z = R*h*(rand(R^2,2)-1/2);
a{8} = Z*[1; 1i];  V = evan(a{8}, 'even');  A{8} = [nq*V, ndqr*V];

for i = 1:length(a)
	% find norms the obvious way
	nrms = pinv(A{i})*Aps;  nrms = sqrt(sum(abs(nrms).^2));  nrms = reshape(nrms, size(x));
	
	% find residuals by taking components in the orthogonal space
	[U,~,~] = svd(A{i});  U = U(:,rank(A{i})+1:end);
	rsdl = U'*Aps;  rsdl = sqrt(sum(abs(rsdl).^2));  rsdl = reshape(rsdl, size(x));
	
	figure
	imagesc(x([1 end]), y([1 end]), rsdl),  colormap gray, axis equal
	hold on,  plot(a{i}, 'ow')
	[C,H] = contour(x,y,nrms, [1.5 5 20 100 300],'-r');  clabel(C,H,'Color','red')
	xlabel 'Re \beta', ylabel 'Im \beta'
	ax = gca;  ax.CLim = [0 1];
	saveTightFigure(sprintf('cohrep%d.pdf', i))
end