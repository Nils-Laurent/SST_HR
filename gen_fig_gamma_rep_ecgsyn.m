close all;

addpath('./paper_code/');

load("data_ecg_syn_N.mat", 'ecgsyn_N', 'Fs', 'hrmean', 'hrstd');

[Lx, N_rep] = size(ecgsyn_N);
prec_bpm = 0.2667; % frequency bin per bpm
max_f = 30;

T_x = (0:(Lx - 1))/Fs;

%% change parameters
GSigs = 1:0.2:3;
N_g = length(GSigs);

std_vec_SST2 = zeros(1, N_g);
mean_vec_SST2 = zeros(1, N_g);
std_vec_SST1 = zeros(1, N_g);
mean_vec_SST1 = zeros(1, N_g);
std_vec_RP = zeros(1, N_g);
mean_vec_RP = zeros(1, N_g);

for nr=1:N_rep
    s_syn = ecgsyn_N(:, nr);
    [X_A_SST2, X_A_SST1, X_A_STFT, N_hsz, BPM_X, Nfft, sigma_w] =...
        ECG_TF(s_syn, Fs, max_f, prec_bpm);

    L_hsz = size(X_A_SST2, 2);

    [W_STFT, W_SST, BPM_comp] = ECG_dictionnary(Fs, Nfft, sigma_w, max_f);
        
%     [RP_hr, RP_mean, RP_std, ~, ~] = hr_from_rr(s_syn, Fs, T_x);
    
%     RP_mean = RP_mean*60 - 29;
%     RP_std = RP_std*60;
    
%     std_vec_RP(:) = std_vec_RP(:) + RP_std;
%     mean_vec_RP(:) = mean_vec_RP(:) + RP_mean;

    for ng=1:N_g
        fprintf("rep = %u/%u, gamma %u/%u \n", nr, N_rep, ng, N_g);
        gSig = GSigs(ng);
        [EMD_T2, ke_T2, LB_T, HB_T, Delta_T] =...
            EMD_ECG_fast(X_A_SST2, W_SST, gSig);
        [EMD_T1, ke_T1, LB_T, HB_T, Delta_T] =...
            EMD_ECG_fast(X_A_SST1, W_SST, gSig);

%         EMDsc_Ismall(N_hsz/Fs, 1:size(EMD_T1, 1), EMD_T1);
%         hold on;
%         plot(N_hsz/Fs, ke_T1, 'c-', 'DisplayName', 'RP HR');
%         plot(N_hsz/Fs, mean(ke_T1)*ones(size(ke_T1)), 'c--', 'DisplayName', 'RP HR');
%         plot(N_hsz/Fs, RP_hr(N_hsz)*60 - 29, 'r-', 'DisplayName', 'RP HR');
%         plot(N_hsz/Fs, RP_mean*ones(size(ke_T1)), 'r--', 'DisplayName', 'RP HR');
%         hold off;
%         
%         EMDsc_Ismall(N_hsz/Fs, 1:size(EMD_T1, 1), EMD_T2);
%         hold on;
%         plot(N_hsz/Fs, ke_T2, 'c-', 'DisplayName', 'RP HR');
%         plot(N_hsz/Fs, mean(ke_T2)*ones(size(ke_T2)), 'c--', 'DisplayName', 'RP HR');
%         plot(N_hsz/Fs, RP_hr(N_hsz)*60 - 29, 'r-', 'DisplayName', 'RP HR');
%         plot(N_hsz/Fs, RP_mean*ones(size(ke_T2)), 'r--', 'DisplayName', 'RP HR');
%         hold off;
        
%         drawnow
%         pause;

        std_vec_SST2(ng) = std_vec_SST2(ng) + std(ke_T2);
        mean_vec_SST2(ng) = mean_vec_SST2(ng) + mean(ke_T2);
        std_vec_SST1(ng) = std_vec_SST1(ng) + std(ke_T1);
        mean_vec_SST1(ng) = mean_vec_SST1(ng) + mean(ke_T1);
    end
end


std_vec_SST2 = std_vec_SST2/N_rep;
mean_vec_SST2 = mean_vec_SST2/N_rep;
std_vec_SST1 = std_vec_SST1/N_rep;
mean_vec_SST1 = mean_vec_SST1/N_rep;
% std_vec_RP = std_vec_RP/N_rep;
% mean_vec_RP = mean_vec_RP/N_rep;

% save("data_fig_gamma_ecgsyn_N.mat",...
%     'N_hsz', 'BPM_X', 'Nfft', 'sigma_w',...
%     'std_vec_SST1', 'mean_vec_SST1', 'std_vec_SST2', 'mean_vec_SST2', 'std_vec_RP', 'mean_vec_RP');

% load("data_fig_gamma_ecgsyn_N.mat");

% [hr, hrm, rpeaks_y, rpeaks_ind] = hr_from_rr(ecgsig, Fs, T_x);
% 
% figure;
% hold on;
% plot(60*hr);
% plot(60*hrm*ones(size(hr)));
% hold off;


%% 2 figures
axisFSZ = 22;
labelSZ = 36;
lenSZ = 700;

figPlot_Ismall("$\gamma$", "components");
legend_Ismall('southwest');
hold on;
plot(GSigs, hrstd*ones(1, N_g), 'DisplayName', '$\sigma_{syn}$',...
    'LineWidth', 2);
plot(GSigs, std_vec_SST1, '--', 'DisplayName', 'std($\mathbf{\hat i}$) FSST',...
    'LineWidth', 2);
plot(GSigs, std_vec_SST2, '-.', 'DisplayName', 'std($\mathbf{\hat i}$) FSST2',...
    'LineWidth', 2);
hold off;
ylim([1.5, max(std_vec_SST2) + 0.5]);
saveas(gcf, "fig_std_ecgsyn", 'epsc');

figPlot_Ismall("$\gamma$", "component index");
legend_Ismall('southeast');
hold on;
% in our case, 60 bpm is the index 31 of the dictionary
plot(GSigs, 31*ones(1, N_g), 'DisplayName', 'HR$_{syn}$',...
    'LineWidth', 2);
plot(GSigs, mean_vec_SST1, '--', 'DisplayName', 'mean($\mathbf{\hat i}$) FSST',...
    'LineWidth', 2);
plot(GSigs, mean_vec_SST2, '-.', 'DisplayName', 'mean($\mathbf{\hat i}$) FSST2',...
    'LineWidth', 2);
hold off;

saveas(gcf, "fig_mean_ecgsyn", 'epsc');