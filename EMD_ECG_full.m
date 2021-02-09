function [EMD_all] = EMD_ECG_full(X_A, We)

addpath('./SST_compare-master/comparison');
addpath('./SST_compare-master/FastEMD');


[~, L_hsz] = size(X_A);
[~, K_e] = size(We);

fprintf("EMD all\n");
EMD_all = zeros(K_e, L_hsz);
for ke=1:K_e
    fprintf("%u/%u ", ke, K_e);
    EMD_all(ke, :) = EMDMatGen(X_A, We(:, ke));
end
fprintf("\n");
end

