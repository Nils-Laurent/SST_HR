close all;

ecg_name = 'ECGt_raw';
signal_ecg = load("test2.mat", ecg_name);
signal_ecg = signal_ecg.(genvarname(ecg_name));

Fs = 1000;
Lx = min(Fs*60, length(signal_ecg));
s_ecg = signal_ecg(1:Lx);

T_x = (0:(Lx-1))/Fs;

% A priori knowledge
bpm_L = 60;
bpm_H = 85;

[N_hat, std_hat, EMD_all] = ECG_emp_SST(s_ecg, Fs, bpm_L, bpm_H);

% N_hat = [30];
% std_hat = [3.63487355824790];