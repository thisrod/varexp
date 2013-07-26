
beginfig(1)

z0 = (0,0);  z1 = 30mm*up;  z2 = 50mm*right;
z10 = 1/4[z0,z1];  y11 = y10; x11 = x2;
z20 = 3/4[z0,z1];  y21 = y20; x21 = x2;
z16 = (10mm,0);  x17 = x16; y17 = y1;
path locus;
z31 = z20+(7mm, -3mm);
x32 = x16;  y32 = y10;
z33 = z11+(-5mm, -3mm);
locus = z20{right}..{z32-z31}z31--z32--z33;

drawdblarrow z1--z0--z2;
draw z10--z11 dashed withdots;
draw z20--z21 dashed withdots;
label.lft("residual", 1/2[z0,z1]);
label.bot(btex solution norm etex, 1/2[z0,z2]);
draw z10--z11 dashed withdots;
label.rt(btex $\|e\|$ etex, z11);
draw z20--z21 dashed withdots;
label.rt(btex norm of rhs etex, z21);
draw z16--z17 dashed withdots;
label.top(btex $\|x_0\|$ etex, z17);
draw locus;

endfig;

end