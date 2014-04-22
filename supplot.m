function h = supplot( z )
%supplot Plot a superposition of coherent states
%   supplot(z) plots the superposition z on the argand plane, with weights
%   represented by saturation and phase by hue.

R = 0.5*numel(z);  f = z(1:R);  a = z(R+1:2*R);
ds = exp(f+0.5*abs(a).^2);
dm = max(abs(ds));
for i=1:numel(a)
    h = plot(real(a(i)), imag(a(i)), 'ok', ...
        'MarkerFaceColor', ...
        hsv2rgb([0.5+angle(ds(i))/2/pi abs(ds(i))/dm 1]));
    hold on
end

if nargout==0; clear h; end

end

