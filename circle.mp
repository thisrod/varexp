beginfig(1)

scale = 10mm;
N = 5;

z0 = (0.5,1.5)*scale;

pair u[], aprx;
aprx = (0,0);
for i = 0 upto N-1:
	u[i] = up rotated (i/N*360);
	drawarrow (origin--u[i]) scaled scale;
	aprx := aprx + 1/4*(2*3.141592/N)*(z0 dotprod u[i])*u[i];
endfor

dotlabel.top("exact", z0);
dotlabel.bot("approx", aprx);

endfig;

end