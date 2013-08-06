beginfig(1)

pair u[], v[], w, a[];  w = (-4,3)/5;
l = 10mm;
transform plot;  plot := identity scaled l;
path p;

for h = 3, 10, 30, 40:
	u[h] = right;
	v[h] = up rotated -2h;
	a[h] = whatever*u[h];  w = a[h]+whatever*v[h];
	b := max(1.2, 1.2*length(a[h]));
	z[h] = unitvector(1/2[u[h],v[h]]);
	fill (1.2u[h]..1.2z[h]..1.2v[h]--(-b*v[h])..(-b*z[h])..(-b*u[h])--cycle) 
		transformed plot withcolor 0.3[white,black];
	drawdblarrow ((v[h]--(0,0)--u[h]) scaled 3/2 transformed plot);
	drawarrow ((0,0)--w) transformed plot;
	draw ((0,0)--a[h]) transformed plot withcolor 1/3[blue,black];
	draw (a[h]--w) transformed plot withcolor 1/3[blue,black];
	plot := plot shifted (30mm*down);
endfor

label.bot(btex $|\alpha\rangle$ etex, 3/2l*u[3]);
label.lft(btex $|\beta\rangle$ etex, 3/2l*v[3]);
label.ulft(btex $|\psi\rangle$ etex, l*w);

endfig;

end