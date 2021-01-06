function [stat,pval] = portmanteu(res,h,lags)
%Portmanteu test of no autocorrelation up to lag h
%   Detailed explanation goes here
% save dimensions
[T,K] = size(res);                                                          
% compute covariance matrix
C0 = res'*res/T;
% compute test statistic
stat = 0;
for j=1:h 
    % autocovariance matrix of order j
    Cj = rmmissing([lagmatrix(res,j) res]);
    Cj = Cj(:,1:K)'*Cj(:,K+1:end)/T;   
    stat = stat + trace(Cj'/C0*Cj/C0);
end
stat = T*stat;
pval = chi2cdf(stat,K^2*(h-lags),'upper');
end