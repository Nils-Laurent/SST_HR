function [T_hsz, BPM_X, BPM_comp, ...
    N_hat, std_hat, N_hat2, std_hat2, EMD_all_SST, EMD_all_STFT] =...
    ECG_emp_both(ecg_in, Fs, bpm_L, bpm_H)
s_ecg = ecg_in - mean(ecg_in);
s_ecg = hilbert(s_ecg);

bin_pu_bpm = 0.2667; % frequency bin per bpm
nBpm = Fs/2*60;
max_f = 30;
Nfft = ceil(nBpm*bin_pu_bpm);
k_max = floor(max_f*Nfft/Fs) + 1;

prec_w = 10^(-3);
[g, Lh, sigma_w] = create_gaussian_window_len(Fs, Nfft, 4*Fs - 1, prec_w);
Lw = 2*Lh + 1;

F_vec = (0:(Nfft-1))*Fs/Nfft;

fprintf("Nfft = %u\n", Nfft);
BPM_vec = F_vec*60;

fprintf("SST\n");
gamma = 10^(-6);
hsz = 32; % Hop size (w-overlap), shift
[STFT,~,SST2, n_down] = sst2_down_gauss(s_ecg, sigma_w, Fs, Nfft, hsz, prec_w);
T_hsz = n_down/Fs;

X_A_SST = abs(SST2(1:k_max, :));
X_A_STFT = abs(STFT(1:k_max, :));
F_X = F_vec(1:k_max);
BPM_X = BPM_vec(1:k_max);

% figure;
% imagesc(T_hsz, BPM_X, X_A_SST);
% axis xy;
% return;

[NX, Lx2] = size(X_A_SST); %Number of frequency and time frames

%% dictionary parameters

fcmin = 30/60; % 0.5 (30 bpm)
fcmax = 180/60; % 3 (180 bpm)

param_We.minFF = fcmin;
param_We.deltaFF = (fcmax - fcmin);
param_We.C_HMax = 16;
param_We.prec_bpm = 1;
param_We.max_f = max_f;

%% set SST, STFT dictionary matrices

Env = ones(1, Nfft); % constant envelope
[W_SST, K_e, FF_vec] = set_We_sst(param_We, Fs, Nfft, Env);
[W_STFT, ~, ~] = set_We_stft(param_We, Fs, Nfft, Env, sigma_w);

slope = (fcmax-fcmin)/(K_e-1); 
F_comp = (fcmin + slope*(0:K_e-1));
BPM_comp = 60*F_comp;

%% EMD
[N_hat, std_hat, EMD_all_SST] = ECG_emp_analysis(X_A_SST, W_SST, BPM_comp);
[N_hat2, std_hat2, EMD_all_STFT] = ECG_emp_analysis(X_A_STFT, W_STFT, BPM_comp);

end

