function [W_exc, Wxe0, fe0, Wxed] = model_WExcitation_interp(fcmin, slope, Kxe, fe, F, N, m, Coef_interp)

% Function that model the pattern matrix W of the ECG excitation part
% Inputs :
%   fcmin: min cardiac freq
%   slope
%   Kex : number of components
%   fe : freq array from spectrogram
%   F : number of frequency frames
%   N : window size
%   m : Zero-padding factor
% Outputs : 
%   Wxe0 : pattern matrix W
%   f0 : freq array for display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Interpolation
Fnew = F*Coef_interp;
pas_fe = fe(2)-fe(1);
fenew = [0:pas_fe/Coef_interp:max(fe)];

%%For strips widening
% wh = hamming(N);
wh = w_hamming(N);
M = N*m*5;
whSpec = fftshift(abs(fft(wh, M)))*10^-3;
maxE = max(whSpec);
threshold = 0.1*maxE; 
indexS = find(whSpec>threshold);
L = fix(length(indexS)/2);

gamma = 1;
Wxe0 = zeros(Fnew,Kxe);
Wxed = zeros(Fnew,Kxe);
femax = fenew(end);


while (fcmin*gamma <= femax)
    k = 0;
    fk = gamma*(fcmin + slope*k);
    while (fk < femax && k<=Kxe-1)
      [~, ind] = min(abs(fk - fenew));
      Wxe0(ind,k+1) = 1;
      Wxed(ind,k+1) = 1;
      if((ind+L)<Fnew)
          Wxe0(ind+(-L:L),k+1) = conv(Wxe0(ind+(-L:L),k+1), whSpec(indexS(1):indexS(L*2+1)), 'same');      
      else         
          Wxe0(ind,k+1) = 1;
      end
      k = k+1;
      fk = gamma*(fcmin + slope*k);
    end
    gamma = gamma+1;
end
fe0 = (fcmin + slope*(0:Kxe-1));

W_exc = downsample(Wxe0, Coef_interp);

end


