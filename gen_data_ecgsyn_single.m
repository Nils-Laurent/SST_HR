close all;

sfecg = 1000; % ECG sampling frequency
N_HBeats = 450; % approximate number of heart beats
Anoise = 0.0; % Additive uniformly distributed measurement noise
hrstd = 2; % Standard deviation of heart rate (default : 1)
lfhfratio = 0.5; % LF/HF ratio (default : 0.5)
sfint = 1000; % Internal sampling frequency

% Order of extrema: [P Q R S T]
ti = [-70 -15 0 15 100];
ai = [1.2 -5 30 -7.5 0.75];
bi = [0.25 0.1 0.1 0.1 0.4];

addpath('./ecgsyn/');

Fs = sfecg;
Lx = Fs*30;

% hrmean = 60; % Mean heart rate [bpm]
% [s_syn_init, ipeaks] = ecgsyn(sfecg,N_HBeats,Anoise,hrmean,hrstd,lfhfratio,sfint,ti,ai,bi);
% save("data_ecg_syn_init60.mat", 's_syn_init', 'Fs', 'hrmean', 'hrstd');
% 
% hrmean = 70; % Mean heart rate [bpm]
% [s_syn_init, ipeaks] = ecgsyn(sfecg,N_HBeats,Anoise,hrmean,hrstd,lfhfratio,sfint,ti,ai,bi);
% save("data_ecg_syn_init70.mat", 's_syn_init', 'Fs', 'hrmean', 'hrstd');
% 
% hrmean = 80; % Mean heart rate [bpm]
% [s_syn_init, ipeaks] = ecgsyn(sfecg,N_HBeats,Anoise,hrmean,hrstd,lfhfratio,sfint,ti,ai,bi);
% save("data_ecg_syn_init80.mat", 's_syn_init', 'Fs', 'hrmean', 'hrstd');