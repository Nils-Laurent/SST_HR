function [He, We, Env] = Har_Env_NMF(X, param_We, N_it, Fs, Nfft, Env_in)

[NX, L] = size(X);
[We, Ke, FF_vec] = set_We_sst(param_We, Fs, Nfft, Env_in);

costF = 'EucDist';

parameters.costFunc = costF; % 'EucDist', 'KLDiv', 'ISDiv'
parameters.numIter = N_it;
parameters.numComp = Ke;

parameters.fixW = 1;
parameters.initW = We;
parameters.initH = ones(Ke, L);

[~, He] = NMF(X, parameters);

m0_v = median(max(He, [], 1));
Y = abs(He - m0_v);
[~, arg_min] = min(Y(:));
[ke_0, n_0] = ind2sub(size(He), arg_min);
f0 = FF_vec(ke_0);

Env = zeros(NX, 1);

P_max = floor(param_We.max_f/f0);
for p=1:P_max
    kp = floor(p*f0*Nfft/Fs)+1;
    Env(p) = X(kp, n_0);
end

%% figures
% figure;
% imagesc(We);
% axis xy;
% colorbar;
% hold on;
% plot(n_0, ke_0, 'ro', 'MarkerSize', 30);
% hold off;
% title("W(e)");
% 
% figure;
% imagesc(He);
% axis xy;
% colorbar;
% hold on;
% plot(n_0, ke_0, 'ro', 'MarkerSize', 30);
% hold off;
% title("H(e)");
% 
% figure;
% imagesc(X);
% axis xy;
% colorbar;
% hold on;
% 
% P_max = floor(param_We.max_f/f0);
% for p=1:P_max
%     kp = floor(p*f0*Nfft/Fs)+1;
%     plot(n_0, kp, 'ro', 'MarkerSize', 30);
% end
% 
% hold off;
% title("X");

end