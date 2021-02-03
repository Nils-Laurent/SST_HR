function [g, Lh] = create_gaussian_window(Fs, Nfft, sigma_w, prec)

% prec = 10^(-3);

Lh = floor(Fs*sigma_w*sqrt(-log(prec)/pi))+1;
Lw = 2*Lh + 1;

t=(1:Lw)'/Fs - (Lh + 1)/Fs;
g = exp(-(t/sigma_w).^2 * pi);

if 2*Lh + 1 > Nfft
    fprintf("[Warning] 2*Lh+1 > Nfft, simga_w = %f\n", sigma_w);
end

end