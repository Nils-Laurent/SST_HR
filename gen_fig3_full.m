close all;


load("data_fig3_GSig_ecgsyn.mat");

Fs = sfecg;
Lx = min(Fs*30, length(s_syn_init));
s_syn = s_syn_init(1:Lx);

prec_bpm = 0.2667; % frequency bin per bpm
max_f = 30;

%% change parameters
bpm_L = 30;
bpm_H = 180;
[T_hsz, BPM_X, BPM_comp, ~, ~, ~, ~, EMD_all_SST, EMD_all_STFT] =...
    ECG_emp_both(s_syn, Fs, bpm_L, bpm_H);

EMDsc_Ismall(T_hsz, BPM_comp, EMD_all_STFT);
saveas(gcf, "fig3_EMD_full_STFT_ecgsyn", 'epsc');

EMDsc_Ismall(T_hsz, BPM_comp, EMD_all_SST);
saveas(gcf, "fig3_EMD_full_SST_ecgsyn", 'epsc');