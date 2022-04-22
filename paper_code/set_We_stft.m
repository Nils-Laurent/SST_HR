function [We, Ke, FF_vec] = set_We_stft(param_We, Fs, Nfft, Env, sigma)
%% create harmonic dictionary
% Fs : sampling frequency
% Nfft : number of frequency bins
% max_f : maximum frequency in dictionary
% minFF : minimum of Fundamental Frequency (FF)
% deltaFF : maximum difference between FF values during time
% C_HMax : coefficient on maximal harmonic (less regularity but smaller
%          matrix)

max_f = param_We.max_f;
minFF = param_We.minFF;
deltaFF = param_We.deltaFF;
% C_HMax = param_We.C_HMax;
prec_bpm = param_We.prec_bpm;
Ke = ceil((deltaFF)*60/prec_bpm) + 1;

k_max = floor(max_f*Nfft/Fs) + 1;
Pe = ceil(max_f/minFF) - 1; % number of harmonic (cuts)
% Pe = ceil(max_f/(minFF + deltaFF)) - 1; % number of harmonic (no cuts)

C = minFF;
% A = C_HMax*Fs*(Fs/Nfft)*(1/Pe);
% Ke = round(deltaFF*Fs/A) + 1;

A = deltaFF*Fs/(Ke - 1);
xe = (0:Ke-1)/Fs;
We = zeros(k_max, Ke);

FF_vec = (C + A*xe);

eta_vec = ((1:k_max) - 1)*Fs/Nfft;
for n=1:Ke
    FF_ap_vec = zeros(1, Pe);
    for p=1:Pe
        FF_ap_vec(p) = C*p + A*p*xe(n);
    end
    
    k = 0;
    for eta = eta_vec
        VPe = exp(-pi*sigma^2*(eta - FF_ap_vec).^2);
        k = k + 1;
        We(k, n) = sum(VPe);
    end
end

% for p=1:Pe
%     m_k_vec = floor((C*p + A*p*xe)*Nfft/Fs) + 1;
%     for n=1:Ke
%         k = m_k_vec(n);
%         if (k > k_max)
%             break;
%         end
%         We(k, n) = Env(p);
%     end
% end

end