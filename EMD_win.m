function [EMD_base, ke_L, ke_H, r] = EMD_win(X_A, We, N_hat)

X_win = X_A(:, 1:N_hat);
[~, K_e] = size(We);

EMD_base = zeros(K_e, N_hat);
fprintf("Computing EMD base (%u)\n", K_e);
for ke=1:K_e
    EMD_base(ke, :) = EMDMatGen(X_win, We(:, ke));
end

[~, ke_vec_b] = min(EMD_base(:, 1:N_hat), [], 1);

N_vec = zeros(1, K_e);

for ke=ke_vec_b
    N_vec(ke) = N_vec(ke) + 1;
end
N_vec = N_vec/N_hat;

[ke_M] = floor(median(ke_vec_b));

ke_L = ceil(3*ke_M/4);
ke_H = floor(3*ke_M/2);

r = sum(N_vec(ke_L:ke_H));
% r

% EMDsc_Ismall(EMD_base);
% hold on;
% plot(ke_vec_b);
% hold off;

% figure;
% hold on;
% plot(1:K_e, N_vec);
% plot([ke_L, ke_M, ke_H], [N_vec(ke_L), N_vec(ke_M), N_vec(ke_H)], '--');
% hold on;

end

