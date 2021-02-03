function [X_A_SST, X_A_STFT, T_hsz, BPM_X, Nfft, sigma_w] =...
    ECG_TF(ecg_in, Fs, max_f, prec_bpm)
s_ecg = ecg_in - mean(ecg_in);
s_ecg = hilbert(s_ecg);

nBpm = Fs/2*60;
Nfft = ceil(nBpm*prec_bpm);
k_max = floor(max_f*Nfft/Fs) + 1;

prec_w = 10^(-3);
[~, ~, sigma_w] = create_gaussian_window_len(Fs, Nfft, 4*Fs - 1, prec_w);
% Lw = 2*Lh + 1;

F_vec = (0:(Nfft-1))*Fs/Nfft;

fprintf("Nfft = %u\n", Nfft);
BPM_vec = F_vec*60;

fprintf("SST\n");
gamma = 10^(-6);
hsz = 32; % Hop size (w-overlap), shift
[STFT,~,SST2, n_down] =...
    sst2_down_gauss(s_ecg, sigma_w, Fs, Nfft, hsz, prec_w, gamma);
T_hsz = n_down/Fs;

X_A_SST = abs(SST2(1:k_max, :));
X_A_STFT = abs(STFT(1:k_max, :));
% F_X = F_vec(1:k_max);
BPM_X = BPM_vec(1:k_max);

end

