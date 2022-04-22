close all;

addpath('./paper_code/');

%% 70 bpm
load("data_ecg_syn_init70.mat", 's_syn_init', 'Fs', 'hrmean', 'hrstd');

Lx = min(Fs*30, length(s_syn_init));
s_syn = s_syn_init(1:Lx);

prec_bpm = 0.2667; % frequency bin per bpm
max_f = 30;

%% change parameters
[X_A_SST2, X_A_SST1, X_A_STFT, T_hsz, BPM_X, Nfft, sigma_w] =...
    ECG_TF(s_syn, Fs, max_f, prec_bpm);
[W_STFT, W_SST, BPM_comp] = ECG_dictionnary(Fs, Nfft, sigma_w, max_f);

[EMD_all_SST2] = EMD_ECG_full(X_A_SST2, W_SST);
[EMD_all_SST1] = EMD_ECG_full(X_A_SST1, W_SST);
[EMD_all_STFT] = EMD_ECG_full(X_A_STFT, W_STFT);

% save("data_fig_full_EMD_ecgsyn70.mat",...
%     'X_A_SST2', 'X_A_SST1', 'X_A_STFT', 'T_hsz', 'BPM_X', 'Nfft', 'sigma_w',...
%     'W_STFT', 'W_SST', 'BPM_comp',...
%     'EMD_all_SST2', 'EMD_all_SST1', 'EMD_all_STFT');

% load("data_fig_full_EMD_ecgsyn70.mat");

ke_L = 30;
ke_H = 50;
[min_T2, r_T2, N_T2] = EMD_ECG_emp_analysis(EMD_all_SST2, ke_L, ke_H);
[min_T1, r_T1, N_T1] = EMD_ECG_emp_analysis(EMD_all_SST1, ke_L, ke_H);
[min_V, r_V, ~] = EMD_ECG_emp_analysis(EMD_all_STFT, ke_L, ke_H);

% N_T
fprintf("------------ 70 ------------\n");
fprintf("B hat T2 = %f\n", r_T2);
fprintf("B hat T1 = %f\n", r_T1);
fprintf("B hat V = %f\n", r_V);




%% fig 70 bpm
close all;

axisFSZ = 22;
labelSZ = 36;
lenSZ = 700;

TFRsc(T_hsz, BPM_X, X_A_STFT, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);
TFRsc(T_hsz, BPM_X, X_A_SST1, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);
TFRsc(T_hsz, BPM_X, X_A_SST2, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);

TFRsc(1:size(W_STFT, 2), BPM_X, W_STFT,...
    "component index", "bpm", 1, axisFSZ, labelSZ, lenSZ);
TFRsc(1:size(W_SST, 2), BPM_X, W_SST,...
    "component index", "bpm", 1, axisFSZ, labelSZ, lenSZ);

EMDsc_Ismall(T_hsz, 1:size(EMD_all_STFT, 1), EMD_all_STFT);
hold on;
plot(T_hsz, min_V, 'w', 'DisplayName', 'EMD min STFT',...
    'LineWidth', 2);
hold off;
% legend_Ismall('northeast');

EMDsc_Ismall(T_hsz, 1:size(EMD_all_SST1, 1), EMD_all_SST1);
hold on;
plot(T_hsz, min_T1, 'w', 'DisplayName', 'EMD min SST',...
    'LineWidth', 2);
hold off;
% legend_Ismall('northeast');

EMDsc_Ismall(T_hsz, 1:size(EMD_all_SST2, 1), EMD_all_SST2);
hold on;
plot(T_hsz, min_T2, 'w', 'DisplayName', 'EMD min SST',...
    'LineWidth', 2);
hold off;
% legend_Ismall('northeast');
pause;
close all;

%% 80 bpm
load("data_ecg_syn_init80.mat", 's_syn_init', 'Fs', 'hrmean', 'hrstd');

Lx = min(Fs*30, length(s_syn_init));
s_syn = s_syn_init(1:Lx);

prec_bpm = 0.2667; % frequency bin per bpm
max_f = 30;

%% change parameters
[X_A_SST2, X_A_SST1, X_A_STFT, T_hsz, BPM_X, Nfft, sigma_w] =...
    ECG_TF(s_syn, Fs, max_f, prec_bpm);
[W_STFT, W_SST, BPM_comp] = ECG_dictionnary(Fs, Nfft, sigma_w, max_f);

[EMD_all_SST2] = EMD_ECG_full(X_A_SST2, W_SST);
[EMD_all_SST1] = EMD_ECG_full(X_A_SST1, W_SST);
[EMD_all_STFT] = EMD_ECG_full(X_A_STFT, W_STFT);

% save("data_fig_full_EMD_ecgsyn80.mat",...
%     'X_A_SST2', 'X_A_SST1', 'X_A_STFT', 'T_hsz', 'BPM_X', 'Nfft', 'sigma_w',...
%     'W_STFT', 'W_SST', 'BPM_comp',...
%     'EMD_all_SST2', 'EMD_all_SST1', 'EMD_all_STFT');

% load("data_fig_full_EMD_ecgsyn80.mat");

ke_L = 40;
ke_H = 60;
[min_T2, r_T2, N_T2] = EMD_ECG_emp_analysis(EMD_all_SST2, ke_L, ke_H);
[min_T1, r_T1, N_T1] = EMD_ECG_emp_analysis(EMD_all_SST1, ke_L, ke_H);
[min_V, r_V, ~] = EMD_ECG_emp_analysis(EMD_all_STFT, ke_L, ke_H);

% N_T
fprintf("------------ 80 ------------\n");
fprintf("B hat T2 = %f\n", r_T2);
fprintf("B hat T1 = %f\n", r_T1);
fprintf("B hat V = %f\n", r_V);

%% fig 80 bpm
close all;

axisFSZ = 22;
labelSZ = 36;
lenSZ = 700;

TFRsc(T_hsz, BPM_X, X_A_STFT, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);
TFRsc(T_hsz, BPM_X, X_A_SST1, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);
TFRsc(T_hsz, BPM_X, X_A_SST2, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);

TFRsc(1:size(W_STFT, 2), BPM_X, W_STFT,...
    "component index", "bpm", 1, axisFSZ, labelSZ, lenSZ);
TFRsc(1:size(W_SST, 2), BPM_X, W_SST,...
    "component index", "bpm", 1, axisFSZ, labelSZ, lenSZ);

EMDsc_Ismall(T_hsz, 1:size(EMD_all_STFT, 1), EMD_all_STFT);
hold on;
plot(T_hsz, min_V, 'w', 'DisplayName', 'EMD min STFT',...
    'LineWidth', 2);
hold off;
% legend_Ismall('northeast');

EMDsc_Ismall(T_hsz, 1:size(EMD_all_SST1, 1), EMD_all_SST1);
hold on;
plot(T_hsz, min_T1, 'w', 'DisplayName', 'EMD min SST',...
    'LineWidth', 2);
hold off;
% legend_Ismall('northeast');

EMDsc_Ismall(T_hsz, 1:size(EMD_all_SST2, 1), EMD_all_SST2);
hold on;
plot(T_hsz, min_T2, 'w', 'DisplayName', 'EMD min SST',...
    'LineWidth', 2);
hold off;
% legend_Ismall('northeast');
