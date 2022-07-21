function Dist = EMDMatGen(X,W)
[Nfft, L] = size(X);
Dist = zeros(1, L);
for n = 1:L
    dX = X(:,n);
    sX = sum(dX);
    if sX > 0
        dX = dX/sX; % normalized 
    else
        dX(1) = 1;
    end
%     dW = zeros(Nfft,1);
%     dW(round(Nfft*W(n)/fmax)) = 1;
    dW = W/sum(W);
    Dist(n) = distVec(dX,dW);
end