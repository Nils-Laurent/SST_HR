close all;

ecg_name = 'ECGt_raw';
signal_ecg = load("test2.mat", ecg_name);
signal_ecg = signal_ecg.(genvarname(ecg_name));

Fs = 1000;
Lx = min(Fs*60, length(signal_ecg));
s_ecg = signal_ecg(1:Lx);

T_x = (0:(Lx-1))/Fs;

%% EMD
% [T_hsz, BPM_X, BPM_comp, N_hat, std_hat, N_hat2, std_hat2, EMD_all_SST, EMD_all_STFT] =...
%     ECG_emp_both(s_ecg, Fs);
% save("data_fig2_EMD_full_ecgreal.mat", 'T_hsz', 'BPM_X', 'BPM_comp',...
%     'N_hat', 'std_hat', 'N_hat2', 'std_hat2', 'EMD_all_SST', 'EMD_all_STFT');


load('data_fig2_EMD_full_ecgreal.mat');

EMDsc_Ismall(T_hsz, BPM_comp, EMD_all_STFT);
saveas(gcf, "fig2_EMD_full_STFT_ecgreal", 'epsc');

EMDsc_Ismall(T_hsz, BPM_comp, EMD_all_SST);
saveas(gcf, "fig2_EMD_full_SST_ecgreal", 'epsc');

