function [W] = w_hamming(Lw)
 a0 = 0.53836;
 Lh = ceil(Lw/2);
 t = (-Lh:Lh)/Lw;
 W = a0 + (1 - a0)*cos(2*pi*t);
end

