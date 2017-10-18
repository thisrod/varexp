% stability of the exact semi-implicit formula

beginfig(1)

u := 4pc;

label(btex $\circ$ etex, (0,0));
z1 = u*dir 220;  dotlabel.llft(btex $s_n$ etex, z1);
z2 = u*dir 100;  dotlabel.top(btex $s_{n+1}$ etex, z2);
path p;
p := fullcircle cutafter ((0,0)--dir 220) cutbefore ((0,0)--dir 40);
p := reverse p scaled 2u;
z3 = point 0 of p;
z4 = z3 + (arclength p) * unitvector(direction 0 of p);
z5 = point infinity of p;
label.urt(btex $\exp(-ihH)s_n$ etex, z5);
drawarrow z3--z4;
label.llft(btex $-ihHs_n$ etex, 1/2[z3,z4]);
draw fullcircle scaled (2u) withcolor 0.6 white;
drawarrow p;
drawarrow z1--z2;
picture f;
f = thelabel(btex $hf_{n+{1\over2}}$ etex, 1/2[z1,z2]);
unfill bbox f;
draw f;

endfig

end