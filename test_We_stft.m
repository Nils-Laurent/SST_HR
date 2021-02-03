close all;

Fs = 1000;

bin_pu_bpm = 0.2667; % frequency bin per bpm
nBpm = Fs/2*60;
max_f = 40;
Nfft = ceil(nBpm*bin_pu_bpm);
k_max = floor(max_f*Nfft/Fs) + 1;

sigma_w = 1.34741677071870;

fcmin = 30/60; % 0.5 (30 bpm)
fcmax = 180/60; % 3 (180 bpm)

param_We.minFF = fcmin;
param_We.deltaFF = (fcmax - fcmin);
param_We.C_HMax = 16;
param_We.prec_bpm = 1;
param_We.max_f = max_f;

fprintf("create dictionary\n");

Env = ones(1, Nfft); % constant envelope
[W_STFT, K_e, FF_vec] = set_We_stft(param_We, Fs, Nfft, Env, sigma_w);

slope = (fcmax-fcmin)/(K_e-1); 
F_comp = (fcmin + slope*(0:K_e-1));
BPM_comp = 60*F_comp;

% close all;
figure;
imagesc(BPM_comp, 1:end, W_STFT);
axis xy;
colorbar;
title("W (STFT)");