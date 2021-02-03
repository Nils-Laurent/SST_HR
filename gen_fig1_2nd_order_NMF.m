close all;

%% global parameters
Fs= 64;

Lx = 4096;
Tx = (0:Lx - 1)/Fs;

%% signal def

P = 4;
A = 1;
a = 1/16;
b = 4;
c = 0;

sig_f = zeros(1, Lx);
for p=1:P
    sig_f = sig_f + A*exp(2*1i*pi*p*(a*Tx.^2/2 + b*Tx + c));
end

sigma_f = 1/sqrt(a);

%% TFR
Nfft = 1024;

prec = 10^(-3);
[g, Lh] = create_gaussian_window(Fs, Nfft, sigma_f, prec);
Lw = 2*Lh + 1;

% max_f = 30;
% k_max = floor(max_f*Nfft/Fs) + 1;
F_vec = (0:(Nfft-1))*Fs/Nfft;
k_max = 512;
max_f = (k_max-1)*Fs/Nfft;


fprintf("SST\n");
gamma = 10^(-6);
hsz = 4; % Hop size (w-overlap), shift
[STFT, ~, SST2, n_down] = sst2_down_gauss(sig_f, sigma_f, Fs, Nfft, hsz, prec);
T_hsz = n_down/Fs;

X_A_SST = abs(SST2(1:k_max, :));
X_A_STFT = abs(STFT(1:k_max, :));
F_X = F_vec(1:k_max);

NX = size(X_A_SST, 1); %Number of frequency frames
Lx2 = size(X_A_SST, 2); %Number of time frames

figure;
imagesc(T_hsz, F_X, X_A_SST);
axis xy;
colorbar;
title("X (abs SST)");

% figure;
% imagesc(T_hsz, F_X, X_A_STFT);
% axis xy;
% colorbar;
% title("X (abs STFT)");
% pause;
% close all;

%% set SST dictionary matrix


fcmin = 4;
fcmax = 8;

param_We.minFF = fcmin;
param_We.deltaFF = (fcmax - fcmin);
param_We.C_HMax = 1;
param_We.max_f = max_f;

fprintf("model_WExcitation_interp\n");

Env = ones(1, Nfft); % constant envelope
[W_SST, K_e, FF_vec] = set_We_sst(param_We, Fs, Nfft, Env);

slope = (fcmax-fcmin)/(K_e-1); 
F_comp = (fcmin + slope*(0:K_e-1));

%% set STFT dictionary matrix

m = 8; %Zero-padding factor

% K_e = 151;

Coef_interp = 10;

[W_STFT, ~, ~, ~] = ...
    model_WExcitation_interp(fcmin, slope, K_e, F_X, NX, Lw, m, Coef_interp);

figure;
imagesc(F_comp, F_X, W_SST);
axis xy;
colorbar;
title("W (SST)");

% figure;
% imagesc(F_comp, F_X, W_STFT);
% axis xy;
% colorbar;
% title("W (STFT)");
% pause;
% close all;

%% minimize (NMF)

costF = 'EucDist';
nbIt = 100;

parameters.costFunc = costF; % 'EucDist', 'KLDiv', 'ISDiv'
parameters.numIter = nbIt;
parameters.numComp = K_e;

parameters.fixW = 1;
parameters.initH = ones(K_e, Lx2);

% ------------ SST
parameters.initW = W_SST;
A = eps+sum(X_A_SST(:));
[~, H_SST] = NMF(X_A_SST, parameters);

V_SST = A*(W_SST*H_SST);

% ------------ STFT
parameters.initW = W_STFT;
A = eps+sum(X_A_STFT(:));
[~, H_STFT] = NMF(X_A_STFT, parameters);

V_STFT = A*(W_STFT*H_STFT);

% close all;

% figure;
% imagesc(T_hsz, F_X, V_STFT);
% axis xy;
% colorbar;
% title("W*H STFT");
% 
% figure;
% imagesc(T_hsz, F_comp, H_STFT);
% axis xy;
% colorbar;
% title("H STFT");

figure;
imagesc(T_hsz, F_X, V_SST);
axis xy;
colorbar;
title("W*H SST");

figure;
imagesc(T_hsz, F_comp, H_SST);
axis xy;
colorbar;
title("H SST");