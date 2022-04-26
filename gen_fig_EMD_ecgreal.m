error("real ECG signal is not provided");

close all;

addpath('./paper_code/');

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

[X_A_SST2, X_A_SST1, X_A_STFT, T_hsz, BPM_X, Nfft, sigma_w] =...
    ECG_TF(s_ecg, Fs, max_f, prec_bpm);
[W_STFT, W_SST, BPM_comp] = ECG_dictionnary(Fs, Nfft, sigma_w, max_f);
[EMD_T, ke_T, LB_T, HB_T] = EMD_ECG_fast(X_A_SST1, W_SST, gSig);

% save("data_fig_EMD_ecgreal.mat",...
%     'T_hsz', 'BPM_X', 'BPM_comp', 'X_A_STFT', 'X_A_SST',...
%     'EMD_T', 'ke_T', 'LB_T', 'HB_T');
% load('data_fig_EMD_ecgreal.mat');

std_vec_SST = std(ke_T);

%% 1 figure

EMDsc_Ismall(T_hsz, 1:size(EMD_T, 1), EMD_T);
plotEMDmin_Ismall(T_hsz, ke_T, 'w', 'HR detection');
% legend_Ismall();
saveas(gcf, "fig_EMD_SST_ecgreal", 'epsc');

