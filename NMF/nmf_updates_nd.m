function [H_excitation, W_filter, H_filter, J] = nmf_updates_nd(X, W_excitation, H_excitation_init, ...
                                                                                        W_filter_init, H_filter_init, ... 
                                                                                        t_spectro, f_spectro, ...
                                                                                        cardiac_freq_array, ...
                                                                                        gammaSmooth, nbIt, flagPlot, num_fig)

                                                                                    
% Function to perform NMF with multiplicative update and smoothness penalization on W_filter
% beta = 2 => Euclidian distance

Kf = size(W_filter_init, 2);
Ke = size(W_excitation, 2);
N = size(H_filter_init, 2);

Val = (gammaSmooth/(N*Ke*Kf));
Val2 = 1/(N*Ke*Kf);

V_excitationtild = W_excitation*H_excitation_init;
V_filtertild = W_filter_init*H_filter_init;
Xtild = V_excitationtild.*V_filtertild;

JData(1) = computeBetaDivergence(X, Xtild);
JS(1) = computeSmoothness(H_excitation_init, W_filter_init, H_filter_init, Val);
J(1) = JData(1) + JS(1);

if flagPlot == 1
        figure(num_fig); clf
        subplot(221)
            semilogy(1, JData(1), 'b+'); ylabel('Loss function of data and smoothness penalization'); hold on
            semilogy(1, JS(1), 'r+'); hold on
            grid on
        subplot(223)
       
        subplot(222)
            semilogy(1, J(1), 'm+'); ylabel('LF of constrained data'); hold on
            grid on
        subplot(224)
        
        drawnow
end

indexRowsZeroWe = find(prod(W_excitation==0,2));

We = W_excitation;
for indIt = 2:nbIt    

    % Update He
    [H_excitation] = updateHe(X, Xtild, ...
                                     W_excitation, H_excitation_init, W_filter_init, H_filter_init, ...
                                     gammaSmooth, Val2, Ke);

    V_excitationtild = We*H_excitation;
    Xtild = V_excitationtild.*V_filtertild;
    

    % Update Hf 
    [H_filter] = updateHf(X, Xtild, ...
                               V_excitationtild, H_excitation, H_filter_init, W_filter_init, gammaSmooth, Val2, Kf);
    
    V_filtertild = W_filter_init*H_filter;
    Xtild = V_excitationtild.*V_filtertild;
   
    
    % Update Wf
    [W_filter] = updateWf(X, Xtild, ...
                               V_excitationtild, indexRowsZeroWe, H_excitation, H_filter, W_filter_init, gammaSmooth, Val2);
    
    V_filtertild = W_filter*H_filter;
    Xtild = V_excitationtild.*V_filtertild;
    
    % Computing cost function
    JData(indIt) = computeBetaDivergence(X, Xtild);
    JS(indIt) = computeSmoothness(H_excitation, W_filter, H_filter, Val);
    J(indIt) = JData(indIt) + JS(indIt);
    
%     figure(10503); clf;
%     for k = 1:Kf
%         subplot(Kf, 2, k); plot(W_filter(:,k)); hold on; plot(W_filter_init(:,k));
%         subplot(Kf, 2, k+Kf); plot(H_filter(k,:)); hold on; plot(H_filter_init(k,:));
%     end
%     
%     figure(10504); clf;
%     imagesc(t_spectro, cardiac_freq_array*60, ...
%     (H_excitation_init ./ (ones(size(H_excitation_init,1),1) * sum(H_excitation_init.^2).^(1/2)))); axis xy; colorbar

    
    
 
%     %Normalization
%     coefNorme = sum(H_excitation.^2,1).^(1/2);
%     H_excitation = H_excitation./(ones(size(H_excitation,1),1)*coefNorme);
%     H_filter = H_filter.*(ones(size(H_filter,1),1)*coefNorme);   
    
    
    % Update parameters
    H_excitation_init = H_excitation;
    H_filter_init = H_filter;
    W_filter_init = W_filter;
    

    
if flagPlot == 1
        figure(num_fig)
        subplot(221)
            semilogy(indIt, JData(indIt), 'b+'); ylabel('Loss function of data and smoothness penalization'); hold on
            semilogy(indIt, JS(indIt), 'r+'); hold on
            grid on
        subplot(223)
            plot(indIt, JData(indIt)-JData(indIt-1), 'b*'); ylabel('Difference'); hold on
            grid on
        subplot(222)
            semilogy(indIt, J(indIt), 'm+'); ylabel('LF of constrained data'); hold on
            grid on
        subplot(224)
            plot(indIt, J(indIt)-J(indIt-1), 'm*'); ylabel('Difference'); hold on
            grid on
        
            
            
            
         
        figure(10505); 
        subplot(311); imagesc(t_spectro, f_spectro, X); axis xy; colorbar; ylabel('Original'); ylim([0 50]); 
        subplot(312); imagesc(t_spectro, f_spectro, (W_excitation*H_excitation_init).*(W_filter_init*H_filter_init)); 
        axis xy;  colorbar; ylabel('Rec'); ylim([0 50]);
        subplot(313); imagesc(t_spectro, f_spectro, X - (W_excitation*H_excitation_init).*(W_filter_init*H_filter_init)); 
        axis xy;  colorbar; ylabel('Residual'); ylim([0 50]);
        
        drawnow
     
end


end

%Normalization
    coefNorme = sum(H_excitation.^2,1).^(1/2);
    H_excitation = H_excitation./(ones(size(H_excitation,1),1)*coefNorme);
    H_filter = H_filter.*(ones(size(H_filter,1),1)*coefNorme);   


end

function J = computeBetaDivergence(X, Y)
    J = sum(sum((X-Y).^2));
end

function JS = computeSmoothness(H_excitation, W_filter, H_filter, Val)

    DeltaW = sum((W_filter(2:end, :)-W_filter(1:end-1,:)).^2,1);
    BetaHf = H_filter.^2;
    Cn = DeltaW*BetaHf;
    JS = (Val/2)*sum(sum(H_excitation.^2,1).*Cn,2);    
    
end

function [H_filter_s] = updateHf(X, Xtild, ...
                               V_excitationtild, H_excitation, H_filter, W_filter, gamma, Val, Kf)

    BetaHe = sum(H_excitation.^2,1);
    DeltaW = sum((W_filter(2:end, :)-W_filter(1:end-1,:)).^2,1);

    
    numer  = W_filter.' * (V_excitationtild.*X);
    denom  = W_filter.' * (V_excitationtild.*Xtild)  + (gamma*Val)*diag(DeltaW)*(ones(Kf,1)*BetaHe.*H_filter);
    H_filter_s      = H_filter .* (numer ./ denom);
    % H_filter_s(H_filter == 0) = 0; % est automatiquement vérifée par
    % l'instruction précédente

end

function [W_filter_s] = updateWf(X, Xtild, ...
                               V_excitationtild, indexRowsZeroWe, H_excitation, H_filter, W_filter, gamma, Val)
    
    BetaHf = H_filter.^2;
    BetaHe = sum(H_excitation.^2,1);
    %DeltaH = sum((ones(Kf,1)*BetaHe).*BetaHf,2);
    DeltaH = sum(BetaHe.*sum(BetaHf,1),2);

    numer = zeros(size(W_filter));
    denom = zeros(size(W_filter));
    F = size(W_filter,1);

    
        %First row
        numer(1,:)      = (V_excitationtild(1,:).*X(1,:)) * H_filter.' + ... 
                                   (gamma*Val)*(W_filter(1,:) + W_filter(1+1,:)).*DeltaH.';
        denom(1,:)      = (V_excitationtild(1,:).*Xtild(1,:)) * H_filter.' + ... 
                                    2*(gamma*Val)*W_filter(1,:)*DeltaH;
                
        %Rows 2 to F-1
        numer(2:F-1,:)      = (V_excitationtild(2:F-1,:).*X(2:F-1,:)) * H_filter.' + ... 
                                (gamma*Val)*(2*W_filter(2:F-1,:) + W_filter((2:F-1)-1,:) + ...
                                W_filter((2:F-1)+1,:))*diag(DeltaH);
        denom(2:F-1,:)      = (V_excitationtild(2:F-1,:).*Xtild(2:F-1,:)) * H_filter.' + ...
                                4*(gamma*Val)*W_filter(2:F-1,:)*DeltaH;

        %Last row
        numer(F,:)      = (V_excitationtild(F,:).*X(F,:)) * H_filter.' + ... 
                                (gamma*Val)*(W_filter(F,:) + W_filter(F-1,:)).*DeltaH.';
        denom(F,:)      = (V_excitationtild(F,:).*Xtild(F,:)) * H_filter.' + ...
                                2*(gamma*Val)*W_filter(F,:)*DeltaH;
 
    
    W_filter_s  = W_filter .* (numer ./ denom);
    
%     W_filter_s(W_filter == 0) = 0;
    W_filter_s(indexRowsZeroWe, :) = 0; 
    
end

function [H_excitation_s] = updateHe(X, Xtild, ...
                                     W_excitation, H_excitation, W_filter, H_filter, ...
                                     gamma, Val, Ke)


    DeltaW = sum((W_filter(2:end, :)-W_filter(1:end-1,:)).^2,1);
    BetaHf = H_filter.^2;
    Cn = DeltaW*BetaHf;
    
    V_filtertild = W_filter * H_filter;
    
    numer  = W_excitation.' * (V_filtertild.*X);
    denom  = W_excitation.' * (V_filtertild.*Xtild)  + (gamma*Val)*H_excitation.*(ones(Ke,1)*Cn);
    H_excitation_s      = H_excitation .* (numer ./ denom);
%     H_excitation_s(H_excitation == 0) = 0;
    
    
end
