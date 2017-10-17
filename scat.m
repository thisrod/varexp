% Q distributions and phase for cat states

global N;  N = 60;  brackets
mesh = -5:0.2:5;  phasespace

cs = nq*evan(1.65*[-1 1], 'even');
cp = sum(cs,2);  cm = cs(:,1) - cs(:,2);

figure
zplot(mesh, mesh, Aps'*cp);
axis image, ylim([-3 3]), saveTightFigure va.pdf

figure
zplot(mesh, mesh, abs(Aps'*cp).^2);
axis image, ylim([-3 3]), saveTightFigure vb.pdf

figure
zplot(mesh, mesh, Aps'*cm);
axis image, ylim([-3 3]), saveTightFigure vc.pdf

figure
zplot(mesh, mesh, abs(Aps'*cm).^2);
axis image, ylim([-3 3]), saveTightFigure vd.pdf
