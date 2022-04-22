function [min_vec, p_hat, N_hat] = EMD_ECG_emp_analysis(EMD_all, ke_L, ke_H)

[~, min_vec] = min(EMD_all, [], 1);
L_min = length(min_vec);

p_hat = sum((min_vec >= ke_L).*(min_vec <= ke_H))/L_min; % Bernoulli estimation

fprintf("compute binomal\n");
N_max = 200;
M_max = floor(N_max/2) - 1;
Y_vec = zeros(1, M_max);
for Mi = 1:M_max
    Ni = 2*Mi + 1;
    Y_vec(Mi) = binocdf(Mi, Ni, p_hat);
end

M_hat = find(Y_vec >= 0.0001, 1, 'last');
N_hat = 2*M_hat + 1;

figure;
plot(Y_vec);

return;
% histogram de ke_vec
N_vec = zeros(1, K_e);

for ke=min_vec
    N_vec(ke) = N_vec(ke) + 1;
end
N_vec = N_vec/sum(N_vec);

% [~, ke_peak] = max(N_vec);
[ke_M] = median(min_vec);

[~, I_zero] = find(N_vec == 0);
D_peak = I_zero - ke_M;
ke_H = min(nonzeros(I_zero.*(D_peak > 0)));
ke_L = max(nonzeros(I_zero.*(D_peak < 0)));

if isempty(ke_H)
    ke_H = K_e;
end
if isempty(ke_L)
    ke_L = 1;
end

% figure;
% hold on;
% plot(1:K_e, N_vec);
% plot(ke_L, N_vec(ke_L), 'o');
% plot(ke_peak, N_vec(ke_peak), 'x');
% plot(ke_H, N_vec(ke_H), 'o');
% hold on;

A1 = 0;
for n=1:L_hsz
    c_ke = min_vec(n);
    A1 = A1 + (ke_L <= c_ke && c_ke <= ke_H); % fix : unit of c_ke
end

p_hat = A1/L_hsz; % Bernoulli estimation

fprintf("compute binomal\n");
N_max = 200;
Y_vec = zeros(1, N_max);
for Ni = 1:N_max
    mid = floor(Ni/2);
    Y_vec(Ni) = binocdf(mid, Ni, p_hat);
end

N_hat = find(Y_vec >= 0.0001, 1, 'last');

[~, ke2_vec] = min(EMD_all(ke_L:ke_H, :), [], 1);
std_hat = std(BPM_comp(ke2_vec));
end

