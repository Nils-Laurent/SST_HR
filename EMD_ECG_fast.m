function [EMD, ke_vec, LB, HB, Delta] = EMD_ECG_fast(X_A, We, gSig)
%% EMD partial
% à implémenter : verif min local

[~, L_hsz] = size(X_A);
[~, K_e] = size(We);

Delta = zeros(1, L_hsz);

LB = zeros(1, L_hsz);
HB = zeros(1, L_hsz);

addpath('./SST_compare-master/comparison');
addpath('./SST_compare-master/FastEMD');

N_hat = 30;

[EMD_base, ke_L, ke_H] = EMD_win(X_A, We, N_hat);

LB(1:N_hat) = ke_L;
HB(1:N_hat) = ke_H;

ke_vec = zeros(1, L_hsz);
[~, ke_vec(1:N_hat)] = min(EMD_base(ke_L:ke_H, :), [], 1);
ke_vec(1:N_hat) = ke_vec(1:N_hat) + ke_L - 1;

fprintf("Computing subset EMD\n", K_e);
EMD = inf(K_e, L_hsz);
EMD(:, 1:N_hat) = EMD_base(:, 1:N_hat);

for n=(N_hat+1):L_hsz
    m = n-1;
    if mod(n, 50) == 0
        fprintf("%u/%u, ", n, L_hsz);
    end

    range_n = ceil(gSig*max(1, std(ke_vec(1:m))));
    Delta(n) = range_n;
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

