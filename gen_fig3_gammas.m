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
% [s_syn_init, ipeaks] = ecgsyn(sfecg,N_HBeats,Anoise,hrmean,hrstd,lfhfratio,sfint,ti,ai,bi);

load("data_fig3_GSig_ecgsyn.mat");

Fs = sfecg;
Lx = min(Fs*30, length(s_syn_init));
s_syn = s_syn_init(1:Lx);

prec_bpm = 0.2667; % frequency bin per bpm
max_f = 30;

%% change parameters
GSigs = 1:0.2:3;
N_g = length(GSigs);

% std_vec_STFT = zeros(1, N_g);
% mean_vec_STFT = zeros(1, N_g);
% std_vec_SST = zeros(1, N_g);
% mean_vec_SST = zeros(1, N_g);
% 
% [X_A_SST, X_A_STFT, T_hsz, BPM_X, Nfft, sigma_w] =...
%     ECG_TF(s_syn, Fs, max_f, prec_bpm);
% 
% L_hsz = size(X_A_SST, 2);
% Delta_STFT = zeros(L_hsz, N_g);
% Delta_SST = zeros(L_hsz, N_g);
% 
% [W_STFT, W_SST, BPM_comp] = ECG_dictionnary(Fs, Nfft, sigma_w, max_f);
% 
% for n=1:N_g
%     fprintf("gamma %u/%u\n", n, N_g);
%     gSig = GSigs(n);
%     [EMD_V, ke_V, LB_V, HB_V, Delta_V] = EMD_ECG_fast(X_A_STFT, W_STFT, gSig);
%     [EMD_T, ke_T, LB_T, HB_T, Delta_T] = EMD_ECG_fast(X_A_SST, W_SST, gSig);
% 
%     Delta_STFT(:, n) = Delta_V;
%     Delta_SST(:, n) = Delta_T;
%     std_vec_STFT(n) = std(BPM_comp(ke_V));
%     mean_vec_STFT(n) = mean(BPM_comp(ke_V));
%     std_vec_SST(n) = std(BPM_comp(ke_T));
%     mean_vec_SST(n) = mean(BPM_comp(ke_T));
% end
% 
% save("data_fig3_GSig_ecgsyn.mat", 's_syn_init',...
%     'X_A_SST', 'X_A_STFT', 'T_hsz', 'BPM_X', 'Nfft', 'sigma_w',...
%     'std_vec_STFT', 'std_vec_SST', 'mean_vec_STFT', 'mean_vec_SST',...
%     'Delta_STFT', 'Delta_SST');

% axisFSZ = 22;
% labelSZ = 36;
% lenSZ = 700;
% 
% TFRsc(T_hsz, BPM_X, X_A_STFT, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);
% saveas(gcf, "fig3_STFT_ecgsyn", 'epsc');
% TFRsc(T_hsz, BPM_X, X_A_SST, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);
% saveas(gcf, "fig3_SST_ecgsyn", 'epsc');
% TFRsc(T_hsz, BPM_X, W_STFT, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);
% saveas(gcf, "fig3_W_STFT_ecgsyn", 'epsc');
% TFRsc(T_hsz, BPM_X, W_SST, "time", "bpm", 1, axisFSZ, labelSZ, lenSZ);
% saveas(gcf, "fig3_W_SST_ecgsyn", 'epsc');

figPlot_Ismall("$\gamma$", "bpm");
legend_Ismall();
hold on;
plot(GSigs, std_vec_STFT, 'DisplayName', 'std($\mathbf{\hat i}$) STFT');
plot(GSigs, std_vec_SST, '--', 'DisplayName', 'std($\mathbf{\hat i}$) SST');
hold off;
saveas(gcf, "fig3_std_ecgsyn", 'epsc');

figPlot_Ismall("$\gamma$", "bpm");
legend_Ismall();
hold on;
plot(GSigs, mean_vec_STFT, 'DisplayName', 'mean($\mathbf{\hat i}$) STFT');
plot(GSigs, mean_vec_SST, '--', 'DisplayName', 'mean($\mathbf{\hat i}$) SST');
hold off;
saveas(gcf, "fig3_mean_ecgsyn", 'epsc');