% close all;

ecg_name = 'ECGt_raw';
signal_ecg = load("test2.mat", ecg_name);
signal_ecg = signal_ecg.(genvarname(ecg_name));

Fs = 1000;
Lx = min(Fs*60, length(signal_ecg));
s_ecg = signal_ecg(1:Lx);

T_x = (0:(Lx-1))/Fs;

prec_bpm = 0.2667; % frequency bin per bpm
max_f = 30;
gSig = 3;

%% EMD

% [X_A_SST, X_A_STFT, T_hsz, BPM_X, Nfft, sigma_w] =...
%     ECG_TF(s_ecg, Fs, max_f, prec_bpm);
% [W_STFT, W_SST, BPM_comp] = ECG_dictionnary(Fs, Nfft, sigma_w, max_f);
% [EMD_V, ke_V, LB_V, HB_V] = EMD_ECG_fast(X_A_STFT, W_STFT, gSig);
% [EMD_T, ke_T, LB_T, HB_T] = EMD_ECG_fast(X_A_SST, W_SST, gSig);

% save("data_fig2_ecgreal.mat",...
%     'T_hsz', 'BPM_X', 'BPM_comp', 'X_A_STFT', 'X_A_SST',...
%     'EMD_V', 'ke_V', 'LB_V', 'HB_V', 'EMD_T', 'ke_T', 'LB_T', 'HB_T');
load('data_fig2_ecgreal.mat');

std_vec_STFT = std(ke_V)
std_vec_SST = std(ke_T)

%% 2 figures
EMDsc_Ismall(T_hsz, 1:size(EMD_V, 1), EMD_V);
plotEMDmin_Ismall(T_hsz, ke_V, 'g', 'HR detection');
legend_Ismall();
saveas(gcf, "fig2_EMD_STFT_ecgreal", 'epsc');

EMDsc_Ismall(T_hsz, 1:size(EMD_T, 1), EMD_T);
plotEMDmin_Ismall(T_hsz, ke_T, 'g', 'HR detection');
legend_Ismall();
saveas(gcf, "fig2_EMD_SST_ecgreal", 'epsc');

