function [VC]=fevd_irf(irf,K,h);
% Structural forecast error variance Decomposition
% using estimated impulse response functions for horizon h 
% Input:
% irf: impulse response function
% K: number of variables
% h: horizon 
%
% Output: rows = horizon; first K cols = effect of first shock
%         on all variables
theta = reshape(irf(1,:), [K,K]);
theta=theta';
% Compute mse(1)
theta2=(theta.*theta);
theta_sum=sum(theta2);
VC=zeros(K,K);
for j=1:K
    VC(j,:)=theta2(j,:)./theta_sum;
end;
% Reshape
VC = reshape(VC', [1, K*K]);
% Then compute FEVD for next horizons
% 
Theta_h=theta2;
for i=2:h
    Theta=reshape(irf(i,:), [K,K]);
    Theta=Theta';
    Theta2=(Theta.*Theta);
    Theta_h = Theta_h+Theta2;
    Theta_sum=sum(Theta_h);
    for j=1:K
       VChelp(j,:)=Theta_h(j,:)./Theta_sum;
    end;
    VC = [ VC; reshape(VChelp', [1,K*K] )];
end;

end
