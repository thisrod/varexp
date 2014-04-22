% Set up phase space grid
[x,y] = meshgrid(mesh);  z = x(:)+1i*y(:);
Aps = nq*evan(z,'even');
DAps = ndqr*evan(z,'even');