%% compaire ECG, cardiac signal estimation (harmonics and envelopes)

close all;

sfecg = 64; % ECG sampling frequency
N_HBeats = 60; % approximate number of heart beats
Anoise = 0.00; % Additive uniformly distributed measurement noise
hrmean = 80; % Mean heart rate (per minute)
hrstd = 1; % Standard deviation of heart rate (default : 1)
lfhfratio = 0.5; % LF/HF ratio (default : 0.5)
sfint = 64; % Internal sampling frequency

% Order of extrema: [P Q R S T]
ti = [-70 -15 0 15 100];
ai = [1.2 -5 30 -7.5 0.75];
bi = [0.25 0.1 0.1 0.1 0.4];

[s_syn, ipeaks] = ecgsyn(sfecg,N_HBeats,Anoise,hrmean,hrstd,lfhfratio,sfint,ti,ai,bi);

Lx = min(4096, length(s_syn));
s_syn = s_syn(1:Lx);

s_syn = s_syn - mean(s_syn);

Fs = sfecg;
T_x = (0:(Lx-1))/Fs;

%% TFR
s_syn = hilbert(s_syn);
Nfft = 1024;

prec = 10^(-3);
[g, Lh, sigma_w] = create_gaussian_window_len(Fs, Nfft, 4*Fs - 1, prec);
[g2, Lh2] = create_gaussian_window(Fs, Nfft, sigma_w, prec);
if Lh ~= Lh2
    error("Lh ~= Lh2");
end
Lw = 2*Lh + 1;


max_f = 30;
k_max = floor(max_f*Nfft/Fs) + 1;
F_vec = (0:(Nfft-1))*Fs/Nfft;

BPM_vec = F_vec*60;

fprintf("SST\n");
gamma = 10^(-6);
hsz = 8; % Hop size (w-overlap), shift
[STFT,~,SST2, n_down] = sst2_down_gauss(s_syn, sigma_w, Fs, Nfft, hsz, prec);
T_hsz = n_down/Fs;

X_A_SST = abs(SST2(1:k_max, :));
% X_A_STFT = abs(STFT(1:k_max, :));
F_X = F_vec(1:k_max);
BPM_X = BPM_vec(1:k_max);

NX = size(X_A_SST, 1); %Number of frequency frames
Lx2 = size(X_A_SST, 2); %Number of time frames

figure;
imagesc(T_hsz, BPM_X, X_A_SST);
axis xy;
colorbar;
title("X (abs SST)");

%% We

fcmin = 0.5; %min cardiac frequency : 30bpm (0.5)
fcmax = 3; %max cardiac frequency : 180bpm (3)

param_We.minFF = fcmin;
param_We.deltaFF = (fcmax - fcmin);
param_We.C_HMax = 2;
param_We.max_f = max_f;

fprintf("model_WExcitation_interp\n");

Env = ones(1, Nfft); % constant envelope
[W_excitation, K_e, FF_vec] = set_We_sst(param_We, Fs, Nfft, Env);

slope = (fcmax-fcmin)/(K_e-1); 
cardiac_freq_array = (fcmin + slope*(0:K_e-1));
BPM_K_e = cardiac_freq_array*60;

%% Min ND
K_f = 2; %number of elements of filter part - mod 2

H_filter_init = ones(K_f, Lx2) + diag(5e-1*rand(1, K_f))*rand(K_f, Lx2);
W_filter_init = ones(NX, K_f) + rand(NX, K_f)*diag(5e-1*rand(1, K_f));
H_excitation_init = ones(K_e, Lx2) + 1e-3*rand(K_e, Lx2);

fprintf("nmf_updates_nd\n");
gamma_s = 1e5;
nbIt= 1e2;

flagPlot = 0;
numFig = 1001;
[H_excitation, W_filter, H_filter, J] =...
    nmf_updates_nd(X_A_SST, W_excitation, H_excitation_init, ...
    W_filter_init, H_filter_init, T_hsz, F_X, cardiac_freq_array,...
    gamma_s, nbIt, flagPlot, numFig);

V = (W_excitation*H_excitation).*(W_filter*H_filter);

figure;
imagesc(T_hsz, BPM_X, V);
axis xy;
colorbar;
title("ND : (We*He).*(Wp*Hp)");

figure;
imagesc(T_hsz, BPM_K_e, H_excitation);
axis xy;
colorbar;
title("ND : He");

%% min NMF toolbox
costF = 'EucDist';

parameters.costFunc = costF; % 'EucDist', 'KLDiv', 'ISDiv'
parameters.numIter = nbIt;
parameters.numComp = K_e;

parameters.fixW = 1;
parameters.initW = W_excitation;
parameters.initH = ones(K_e, Lx2);


A = eps+sum(X_A_SST(:));
[W2, H2] = NMF(X_A_SST, parameters);

V2 = A*(W2*H2);

figure;
imagesc(T_hsz, BPM_X, V2);
axis xy;
colorbar;
title("NMF tb : W*H");

figure;
imagesc(T_hsz, BPM_K_e, H2);
axis xy;
colorbar;
title("NMF tb : H");
