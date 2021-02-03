clear all
close all
load test2 % test2 ou test4 ou test6 %% 3 exemples différents

% fréquence d'échantillonnage
Fs= 1000; %1KHz

% choix du signal
signal_ecg=ECGa;

% duree de signal
signal_ecg=signal_ecg(1:Fs*50);


% Spectrogram Parameters
win = 4096; %4096; %Window size for the spectrogram ~4s
m = 2; %Zero-padding factor
h = 32; %Hop size (w-overlap), shift
nfft = win*2; %Number of FFT points (avant win*8)


% calcul et affichage du spectrogramme pour le signal ECGf
signal_ecg = abs(signal_ecg) - mean(abs(signal_ecg));
[s_spectro,f_spectro,t_spectro,p_spectro] = spectrogram(signal_ecg, hamming(win), win-h, nfft,Fs);
X_ecg=abs(s_spectro).^2;
X_ecg = X_ecg/norm(X_ecg, 'fro');

figure;
imagesc(t_spectro, f_spectro, 10*log10(X_ecg));
ylabel('Fetal ECG'); axis xy; colorbar
xlabel('Time [s]', 'Interpreter', 'latex', 'Fontname', 'Times', 'Fontangle', 'italic', 'Fontsize',18);

% si on limite les fréquences à 50Hz pour une meilleure visualisation de la
% bande d'intérêt
fmax = 50;
X_ecg(f_spectro>fmax, :)=[];
f_spectro(f_spectro>fmax)=[];
figure; 
imagesc(t_spectro, f_spectro, 10*log10(X_ecg));
ylabel('ECG'); axis xy; colorbar
xlabel('Time [s]', 'Interpreter', 'latex', 'Fontname', 'Times', 'Fontangle', 'italic', 'Fontsize',18);



%% dictionnaire Wexcitation

fcmin = 0.5; %min cardiac frequency : 30bpm (0.5)
fcmax = 3; %max cardiac frequency : 240bpm (4) pour le foetus ou 180 (3) pour l'adulte
resolution=1; %bpm : 1 ou 0.5

%   Ke : number of components
Ke=(fcmax*60-fcmin*60)/resolution + 1;
% pente
slope = (fcmax-fcmin)/(Ke-1); 

F=size(X_ecg, 1); % number of frequency frames

%Interpolation
% on va suréchantillonner pour créer les lignes puis sous echantilloner à la fin
% nécessaire pour élargir la largeur des lignes du dictionnaire
Coef_interp = 10;   
Fnew = F*Coef_interp;
pas_fe = f_spectro(2)-f_spectro(1);
fenew = [0:pas_fe/Coef_interp:max(f_spectro)];

%%For strips widening
wh = hamming(win);
M = win*m*5;
whSpec = fftshift(abs(fft(wh, M)))*10^-3;
maxE = max(whSpec);
threshold = 0.8*maxE;  % en changeant le seuil, tu peux modifier le nombre de points que tu récupères et donc la largeur des bandes dans le dictionnaire
indexS = find(whSpec>threshold);
L = fix(length(indexS)/2); 
% avec  threshold = 0.01*maxE; => L=15 => qui va peremettre 3 ou 4 points au final après sous-échantillonnage
% avec  threshold = 0.8*maxE; => L=5 => qui va peremettre 1 ou 2 points au final après sous-échantillonnage

gamma = 1;
Wxe0 = zeros(Fnew,Ke);
Wxed = zeros(Fnew,Ke);
femax = fenew(end);


while (fcmin*gamma <= femax)
    k = 0;
    fk = gamma*(fcmin + slope*k);
    while (fk < femax && k<=Ke-1)
      [~, ind] = min(abs(fk - fenew));
      Wxe0(ind,k+1) = 1;
      Wxed(ind,k+1) = 1;
      if((ind+L)<Fnew && (ind-L)>0)
          Wxe0(ind+(-L:L),k+1) = conv(Wxe0(ind+(-L:L),k+1), whSpec(indexS(1):indexS(L*2+1)), 'same');      
      else         
          Wxe0(ind,k+1) = 1;
      end
      k = k+1;
      fk = gamma*(fcmin + slope*k);
    end
    gamma = gamma+1;
end
% card freq array for display
fe0 = (fcmin + slope*(0:Ke-1));

W_exc = downsample(Wxe0, Coef_interp);

figure(100); clf;
imagesc(fe0*60, f_spectro, W_exc); 
axis xy; colorbar;
ylabel({'frequence'; '[s^{-1}]'});
xlabel({'frequence cardiaque (fc)'; '[bpm]'});

figure(100); clf;
imagesc(W_exc); 
axis xy; colorbar;

figure; plot(W_exc(:,1)); 
hold on; plot(W_exc(:,2)); 

figure; plot(W_exc(:,100),'+-'); 


