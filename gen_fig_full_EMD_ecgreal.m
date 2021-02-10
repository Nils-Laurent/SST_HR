close all;

ecg_name = 'ECGt_raw';
signal_ecg = load("test2.mat", ecg_name);
signal_ecg = signal_ecg.(genvarname(ecg_name));

Fs = 1000;
Lx = min(Fs*60, length(signal_ecg));
s_ecg = signal_ecg(1:Lx);

T_x = (0:(Lx-1))/Fs;

%% EMD
% [X_A_SST, X_A_STFT, T_hsz, BPM_X, Nfft, sigma_w] =...
%     ECG_TF(s_ecg, Fs, max_f, prec_bpm);
% [W_STFT, W_SST, BPM_comp] = ECG_dictionnary(Fs, Nfft, sigma_w, max_f);
% 
% [EMD_all_SST] = EMD_ECG_full(X_A_SST, W_SST);
% [EMD_all_STFT] = EMD_ECG_full(X_A_STFT, W_STFT);
% 
% save("data_fig_full_EMD_ecgreal.mat",...
%     'X_A_SST', 'X_A_STFT', 'T_hsz', 'BPM_X', 'Nfft', 'sigma_w',...
%     'W_STFT', 'W_SST', 'BPM_comp',...
%     'EMD_all_SST', 'EMD_all_STFT');

load('data_fig_full_EMD_ecgreal.mat');

ke_L = 25;
ke_H = 60;
[min_T, r_T, N_T] = EMD_ECG_emp_analysis(EMD_all_SST, ke_L, ke_H);
[min_V, r_V, ~] = EMD_ECG_emp_analysis(EMD_all_STFT, ke_L, ke_H);

N_T
r_T

r_V

pause;

%% 2 figs
close all;

EMDsc_Ismall(T_hsz, 1:size(EMD_all_STFT, 1), EMD_all_STFT);
saveas(gcf, "fig_full_EMD_STFT_ecgreal", 'epsc');

EMDsc_Ismall(T_hsz, 1:size(EMD_all_SST, 1), EMD_all_SST);
saveas(gcf, "fig_full_EMD_SST_ecgreal", 'epsc');

