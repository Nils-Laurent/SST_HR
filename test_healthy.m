%% Test healthy
close all;

%% global parameters
Fs = 1000; % 1KHz
hsz = 32; % Hop size (w-overlap), shift

Nfft = 8*Fs;

[g, Lh, sigma_w] = create_gaussian_window_len(Fs, Nfft, 4*Fs - 1);
% [g, Lh] = create_gaussian_window(Fs, Nfft, sigma_w);
Lw = 2*Lh + 1;

L_env = Nfft;
% TH_amp = zeros(1, L_env);
% TH_amp_tmp = zeros(1, L_env);
% for n=1:L_env
%     TH_amp(n) = (1/n);
%     if n > 3
%         TH_amp_tmp(n) = TH_amp(n-3);
%     else
%         TH_amp_tmp(n) = 1/(2^(3-n));
%     end
% end

max_f = 50;

%% signal def
Lx = Fs*60;
t = (0:Lx - 1)/Fs;

signal_ecg = load("HealthyECG.mat");
signal_ecg = signal_ecg.x;

signal_ecg=signal_ecg(Lx:2*Lx);
Tx = (Lx:(2*Lx-1))/Fs;
in_signal = signal_ecg;

% figure;
% plot(in_signal);

k_max = floor(max_f*Nfft/Fs) + 1;
F_vec = (0:(k_max-1))*Fs/Nfft;

%% SST
fprintf("SST\n");
gamma = 10^(-6);
[STFT,~,SST2, n_down] = sst2_down_gauss(in_signal, sigma_w, Fs, Nfft, hsz);
T_vec = n_down/Fs + Lx/Fs;
SST2_ref = abs(SST2(1:k_max, :));
L = size(SST2_ref, 2);

STFT_ref = abs(STFT(1:k_max, :));

figure;
imagesc(Tx, F_vec, STFT_ref);
axis xy;
colorbar;
title("STFT ref");

figure;
imagesc(Tx, F_vec, SST2_ref);
axis xy;
colorbar;
title("SST2 ref");

% pause;

%% settings

param_We.minFF = 1;
param_We.deltaFF = 1;
param_We.C_HMax = 2;
param_We.max_f = max_f;

Env = ones(1, L_env);
[We, Ke, FF_vec] = set_We_sst(param_We, Fs, Nfft, Env);

%% optimisation
fprintf("NMF SST\n");
N_it = 25;

N_all = 1;
for n=1:N_all
    fprintf(" %u", n);
    Env_prec = Env;
    [He, We, Env] = Har_Env_NMF(SST2_ref, param_We, N_it, Fs, Nfft, Env);

%     figure;
%     imagesc(He);
%     axis xy;
%     colorbar;
%     title("H(e)");
end
fprintf("\n", n);


figure;
imagesc(FF_vec, F_vec, We);
axis xy;
colorbar;
title("W(e)");

figure;
imagesc(T_vec, FF_vec, He);
axis xy;
colorbar;
title("H(e)");

% figure;
% plot(Env_prec(1:20));
% title("Env 20");

% figure;
% imagesc(We*He);
% axis xy;
% colorbar;
% title("W(e)*H(e)");
