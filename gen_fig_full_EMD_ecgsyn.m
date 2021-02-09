close all;

load("data_ecg_syn_init.mat", 's_syn_init', 'Fs', 'hrmean', 'hrstd');

Lx = min(Fs*30, length(s_syn_init));
s_syn = s_syn_init(1:Lx);

prec_bpm = 0.2667; % frequency bin per bpm
max_f = 30;

%% change parameters
% [X_A_SST, X_A_STFT, T_hsz, BPM_X, Nfft, sigma_w] =...
%     ECG_TF(s_syn, Fs, max_f, prec_bpm);
% [W_STFT, W_SST, BPM_comp] = ECG_dictionnary(Fs, Nfft, sigma_w, max_f);
% 
% [EMD_all_SST] = EMD_ECG_full(X_A_SST, W_SST);
% [EMD_all_STFT] = EMD_ECG_full(X_A_STFT, W_STFT);
% 
% save("data_fig_full_EMD_ecgsyn.mat",...
%     'X_A_SST', 'X_A_STFT', 'T_hsz', 'BPM_X', 'Nfft', 'sigma_w',...
%     'W_STFT', 'W_SST', 'BPM_comp',...
%     'EMD_all_SST', 'EMD_all_STFT');

load("data_fig_full_EMD_ecgsyn.mat");

ke_L = 20;
ke_H = 40;
[min_T, r_T, N_T] = EMD_ECG_emp_analysis(EMD_all_SST, ke_L, ke_H);
[min_V, r_V, ~] = EMD_ECG_emp_analysis(EMD_all_STFT, ke_L, ke_H);

N_T
r_T

r_V

pause;
% return;

%% 2 fig
close all;

axisFSZ = 22;
labelSZ = 36;
lenSZ = 700;

TFRsc(T_hsz, BPM_X, X_A_STFT, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);
saveas(gcf, "fig_STFT_ecgsyn", 'epsc');
TFRsc(T_hsz, BPM_X, X_A_SST, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);
saveas(gcf, "fig_SST_ecgsyn", 'epsc');

TFRsc(1:size(W_STFT, 2), BPM_X, W_STFT,...
    "component index", "bpm", 1, axisFSZ, labelSZ, lenSZ);
saveas(gcf, "fig_W_STFT_ecgsyn", 'epsc');
TFRsc(1:size(W_SST, 2), BPM_X, W_SST,...
    "component index", "bpm", 1, axisFSZ, labelSZ, lenSZ);
saveas(gcf, "fig_W_SST_ecgsyn", 'epsc');

EMDsc_Ismall(T_hsz, 1:size(EMD_all_STFT, 1), EMD_all_STFT);
hold on;
plot(T_hsz, min_V, 'g', 'DisplayName', 'EMD min STFT',...
    'LineWidth', 2);
hold off;
legend_Ismall('northeast');
saveas(gcf, "fig_full_EMD_STFT_ecgsyn", 'epsc');

EMDsc_Ismall(T_hsz, 1:size(EMD_all_SST, 1), EMD_all_SST);
hold on;
plot(T_hsz, min_T, 'g', 'DisplayName', 'EMD min SST',...
    'LineWidth', 2);
hold off;
legend_Ismall('northeast');
saveas(gcf, "fig_full_EMD_SST_ecgsyn", 'epsc');
