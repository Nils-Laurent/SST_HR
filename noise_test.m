
Lx = 4096;
Nfft = 1024;

x = randn(1, Lx);
x = x - mean(x);

prec_w = 10^(-3);
[g, Lh, sigma_w] = create_gaussian_window_len(4096, Nfft, 256, prec_w);
Lw = 2*Lh + 1;

fprintf("SST\n");
gamma = 10^(-6);
hsz = 1; % Hop size (w-overlap), shift
[STFT,~,SST2, n_down] = sst2_down_gauss(x, sigma_w, 4096, Nfft, hsz, prec_w);


figure;
imagesc(abs(STFT));
axis xy;
colorbar;
title("X");