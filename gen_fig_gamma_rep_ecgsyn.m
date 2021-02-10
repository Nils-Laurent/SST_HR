close all;

load("data_ecg_syn_N.mat", 'ecgsyn_N', 'Fs', 'hrmean', 'hrstd');

[Lx, N_rep] = size(ecgsyn_N);
prec_bpm = 0.2667; % frequency bin per bpm
max_f = 30;

%% change parameters
GSigs = 1:0.2:3;
N_g = length(GSigs);

std_vec_SST = zeros(1, N_g);
mean_vec_SST = zeros(1, N_g);

for nr=1:N_rep
    s_syn = ecgsyn_N(:, nr);
    [X_A_SST, X_A_STFT, T_hsz, BPM_X, Nfft, sigma_w] =...
        ECG_TF(s_syn, Fs, max_f, prec_bpm);

    L_hsz = size(X_A_SST, 2);

    [W_STFT, W_SST, BPM_comp] = ECG_dictionnary(Fs, Nfft, sigma_w, max_f);

    for ng=1:N_g
        fprintf("rep = %u/%u, gamma %u/%u \n", nr, N_rep, ng, N_g);
        gSig = GSigs(ng);
        [EMD_T, ke_T, LB_T, HB_T, Delta_T] = EMD_ECG_fast(X_A_SST, W_SST, gSig);

%         EMDsc_Ismall(EMD_T);
%         hold on;
%         plot(ke_T);
%         plot(mean(ke_T)*ones(size(ke_T)), '--');
%         hold off;
% %         drawnow
% %         pause;

        std_vec_SST(ng) = std_vec_SST(ng) + std(ke_T);
        mean_vec_SST(ng) = mean_vec_SST(ng) + mean(ke_T);
    end
end

std_vec_SST = std_vec_SST/N_rep;
mean_vec_SST = mean_vec_SST/N_rep;

save("data_fig_gamma_ecgsyn_N.mat",...
    'T_hsz', 'BPM_X', 'Nfft', 'sigma_w',...
    'std_vec_SST', 'mean_vec_SST');

load("data_fig_gamma_ecgsyn_N.mat");

%% 2 figures
axisFSZ = 22;
labelSZ = 36;
lenSZ = 700;

figPlot_Ismall("$\gamma$", "components");
legend_Ismall('southwest');
hold on;
plot(GSigs, hrstd*ones(1, N_g), 'DisplayName', '$\sigma_{syn}$');
plot(GSigs, std_vec_SST, '--', 'DisplayName', 'std($\mathbf{\hat i}$) SST');
hold off;
saveas(gcf, "fig_std_ecgsyn", 'epsc');

figPlot_Ismall("$\gamma$", "component index");
legend_Ismall();
hold on;
% in our case, 60 bpm is the index 31 of the dictionary
plot(GSigs, 31*ones(1, N_g), 'DisplayName', 'HR$_{syn}$');
plot(GSigs, mean_vec_SST, '--', 'DisplayName', 'mean($\mathbf{\hat i}$) SST');
hold off;
saveas(gcf, "fig_mean_ecgsyn", 'epsc');