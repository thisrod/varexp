% Plot examples of states for conventions chapter.

global N; N = 24; brackets
mesh = -5:0.1:5;  phasespace
ht = 3.7;

figure		% vacuum state
ax = zplot(mesh, mesh, Aps'*nq*evan(0, 'even'));
ax.XTick = -4:4:4;  ax.YTick = -4:4:4;
ax.Units = 'centimeters';  ax.Position = [0.5 0.5 ht ht];
xlabel 'Re \alpha',  ylabel 'Im \alpha'
axis image
saveTightFigure twoa.pdf

figure
ax = zplot(mesh, mesh, Aps'*nq*evan(3+2i, 'even')); hold on
plot(0,0,'ow')
ax.XTick = -4:2:4;  ax.YTick = -4:2:4;
ax.Units = 'centimeters';  ax.Position = [0.5 0.5 ht ht];
axis off, axis image
saveTightFigure twob.pdf

figure
ax = zplot(mesh, mesh, Aps'*[zeros(5,1); 1; zeros(N-5,1)]); hold on
plot(0,0,'ow')
ax.XTick = -4:2:4;  ax.YTick = -4:2:4;
ax.Units = 'centimeters';  ax.Position = [0.5 0.5 ht ht];
axis off, axis image
saveTightFigure twoc.pdf

r = linspace(0.75,1,10);
theta = linspace(0, 2*pi, 100);
[rg, thg] = meshgrid(r,theta);
[x,y] = pol2cart(thg,rg);
figure, ax = axes;
pcolor(x,y,thg);
colormap(phase(128));
ax.Units = 'centimeters';  ax.Position = [0.5 0.5 ht ht];
axis off, axis image
text(0.85,0,'1'), text(-0.95,0,'-1'), text(0,0.85,'i'), text(0,-0.9,'-i')
shading flat;
saveTightFigure twod.pdf

figure, ax = axes;
image([113:127 0:128 1:15]'), colormap(phase(129))
ax.YDir = 'normal';  ax.XTick = [];
ax.YTick = 16+32*(0:4);  ax.YTickLabel = {'1', 'i', '-1', '-i', '1'};
ax.YAxisLocation = 'right';	% saveTightFigure clips the left side
ax.Units = 'centimeters';  ax.Position = [0.5 0.5 0.6 ht];
