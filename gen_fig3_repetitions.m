close all;

sfecg = 1000; % ECG sampling frequency
N_HBeats = 450; % approximate number of heart beats
Anoise = 0.0; % Additive uniformly distributed measurement noise
hrmean = 60; % Mean heart rate [bpm]
hrstd = 2; % Standard deviation of heart rate (default : 1)
lfhfratio = 0.5; % LF/HF ratio (default : 0.5)
sfint = 1000; % Internal sampling frequency

% Order of extrema: [P Q R S T]
ti = [-70 -15 0 15 100];
ai = [1.2 -5 30 -7.5 0.75];
bi = [0.25 0.1 0.1 0.1 0.4];

addpath('./ecgsyn/');
[s_syn_init, ipeaks] = ecgsyn(sfecg,N_HBeats,Anoise,hrmean,hrstd,lfhfratio,sfint,ti,ai,bi);

Fs = sfecg;
Lx = min(Fs*30, length(s_syn_init));
s_syn = s_syn_init(1:Lx);

T_x = (0:(Lx-1))/Fs;

%% varying noise level
SNRs = inf;
N_snr = length(SNRs);
std_vec_STFT = zeros(1, N_snr);
mean_vec_STFT = zeros(1, N_snr);
std_vec_SST = zeros(1, N_snr);
mean_vec_SST = zeros(1, N_snr);

for n=1:N_snr
    snr = SNRs(n);
    
    noise = randn(1, Lx);
    s_noise = sigmerge(s_syn, noise', snr);
%     [T_hsz, BPM_X, BPM_comp, R_STFT, R_SST] =...
%         ECG_TF_cmp(s_noise, Fs);

    std_vec_STFT(n) = std(R_STFT.CVec);
    mean_vec_STFT(n) = mean(R_STFT.CVec);
    std_vec_SST(n) = std(R_SST.CVec);
    mean_vec_SST(n) = mean(R_SST.CVec);
end

TFRsc_Ismall(T_hsz, BPM_X, R_SST.ASST);
TFRsc_Ismall(T_hsz, BPM_X, R_STFT.ASTFT);
EMDsc_Ismall(T_hsz, BPM_comp, R_SST.EMD);
EMDsc_Ismall(T_hsz, BPM_comp, R_STFT.EMD);
