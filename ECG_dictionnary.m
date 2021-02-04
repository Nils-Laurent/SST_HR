function [W_STFT, W_SST, BPM_comp] = ECG_dictionnary(Fs, Nfft, sigma_w, max_f)


fcmin = 30/60; % 0.5 (30 bpm)
fcmax = 180/60; % 3 (180 bpm)

param_We.minFF = fcmin;
param_We.deltaFF = (fcmax - fcmin);
param_We.C_HMax = 16;
param_We.prec_bpm = 1;
param_We.max_f = max_f;

%% set SST, STFT dictionary matrices

Env = ones(1, Nfft); % constant envelope
[W_SST, K_e, FF_vec] = set_We_sst(param_We, Fs, Nfft, Env);
[W_STFT, ~, ~] = set_We_stft(param_We, Fs, Nfft, Env, sigma_w);

slope = (fcmax-fcmin)/(K_e-1); 
F_comp = (fcmin + slope*(0:K_e-1));
BPM_comp = 60*F_comp;

end

