beginfig(1)

scale = 30mm;
N = 5;
color junk;  junk = 0.5 white;

pair u[], a, b, c;

for i = 0 upto N-1:
	u[i] = up rotated (i/N*180);
	draw (-u[i]--u[i]) scaled scale withcolor junk;
endfor

z1 = (0.7,0.5);  z2 = (0.6,0.6);

for t = z1, z2:
	label("+", t*scale);
	a := (0,0);
	for i = 0 upto N-1:
		b := a;
		a := a + (2/N)*(t dotprod u[i])*u[i];
		draw (b--a) scaled scale;
	endfor
endfor

z11 = whatever*u3 = z1 + whatever*u4;
draw (origin--z11--z1) scaled scale withcolor 1/2[red,white];
z12 = whatever*u3 = z2 + whatever*u4;
draw (origin--z12--z2) scaled scale withcolor 1/2[red,white];


endfig;

end