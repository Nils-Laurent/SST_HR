function [EMD, ke_vec, LB, HB] = EMD_ECG_fast(X_A, We, gSig)
%% EMD partial
% à implémenter : verif min local

[~, L_hsz] = size(X_A);
[~, K_e] = size(We);

LB = zeros(1, L_hsz);
HB = zeros(1, L_hsz);

addpath('./SST_compare-master/comparison');
addpath('./SST_compare-master/FastEMD');

N_hat = 30;
% std_hat = 3.70364723691537; % bpm
% scope = ceil(2*bin_p_bpm*std_hat);

EMD_base = zeros(K_e, N_hat);
fprintf("Computing EMD base (%u)\n", K_e);
for ke=1:K_e
    EMD_base(ke, :) = EMDMatGen(X_A(:, 1:N_hat), We(:, ke));
end

[~, ke_vec_b] = min(EMD_base, [], 1);

N_vec = zeros(1, K_e);

for ke=ke_vec_b
    N_vec(ke) = N_vec(ke) + 1;
end
N_vec = N_vec/sum(N_vec);


[~, ke_peak] = max(N_vec(2:(K_e - 1)));
ke_peak = ke_peak + 1;

[~, I_zero] = find(N_vec == 0);
D_peak = I_zero - ke_peak;
ke_H = min(nonzeros(I_zero.*(D_peak > 0)));
ke_L = max(nonzeros(I_zero.*(D_peak < 0)));

if isempty(ke_H)
    ke_H = K_e;
end
if isempty(ke_L)
    ke_L = 1;
end

LB(1:N_hat) = ke_L;
HB(1:N_hat) = ke_H;

ke_vec = zeros(1, L_hsz);
[~, ke_vec(1:N_hat)] = min(EMD_base(ke_L:ke_H, :), [], 1);
ke_vec(1:N_hat) = ke_vec(1:N_hat) + ke_L - 1;

% figure;
% hold on;
% plot(1:K_e, N_vec);
% plot(ke_L, N_vec(ke_L), 'o');
% plot(ke_peak, N_vec(ke_peak), 'x');
% plot(ke_H, N_vec(ke_H), 'o');
% hold on;

fprintf("Computing subset EMD\n", K_e);
EMD = inf(K_e, L_hsz);
EMD(:, 1:N_hat) = EMD_base(:, 1:N_hat);

for n=(N_hat+1):L_hsz
    m = n-1;
    if mod(n, 50) == 0
        fprintf("%u/%u, ", n, L_hsz);
    end

    range_n = gSig*max(1, ceil(std(ke_vec(1:m))));
    km = ke_vec(m);
    LB(n) = max(1, km - range_n);
    HB(n) = min(K_e, km + range_n);
    for ke=LB(n):HB(n)
        EMD(ke, n) = EMDMatGen(X_A(:, n), We(:, ke));
    end
    [~, v_arg] = min(EMD(LB(n):HB(n), n));
    ke_vec(n) = v_arg + LB(n) -1;
end
fprintf("\n");

end

