close all;

ecg_name = 'ECGt_raw';
signal_ecg = load("test2.mat", ecg_name);
signal_ecg = signal_ecg.(genvarname(ecg_name));

Fs = 1000;
Lx = min(Fs*60, length(signal_ecg));
s_ecg = signal_ecg(1:Lx);

T_x = (0:(Lx-1))/Fs;


%% EMD
% [T_hsz, BPM_X, BPM_comp, R_STFT, R_SST] =...
%     ECG_TF_cmp(s_ecg, Fs);
% save("fig2_data_ecgreal.mat", 'T_hsz', 'BPM_X', 'BPM_comp', 'R_STFT', 'R_SST');


load('data_fig2_ecgreal.mat');

std_vec_STFT = std(R_STFT.CVec)
std_vec_SST = std(R_SST.CVec)

save_fig2_EMD(T_hsz, BPM_X, BPM_comp,...
    R_STFT, R_SST, "ecgreal");

