close all;

%% global parameters
Fs = 1000; % 1KHz

% Choose between : ECGa(_raw), ECGt(_raw), ECGf
ecg_name = 'x';
signal_ecg = load("HealthyECG.mat", ecg_name);
signal_ecg = signal_ecg.(genvarname(ecg_name));

figure;
plot(signal_ecg(1:8000));
return;