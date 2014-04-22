function h = wplot(x, y, w)
%zplot Plot a function of a complex variable as a color image
%   x and y are as for image.  w is a matrix of samples f(x+iy).
%   The modulus of f is plotted on the Argand plane as saturation, 
%   its argument as hue.  Colors might be excessively quantised.

dat = floor(64*abs(w)/max(abs(w(:))));
CV = ind2rgb(dat, gray(64));
dat = floor(32*(1+angle(w)/pi));
CHS = ind2rgb(dat, hsv(64));
h = image(x, y, CV.*(CHS-1)+1);
if nargout==0; clear h; end

end

