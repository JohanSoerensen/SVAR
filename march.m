function [test]=march(res,lags,K)
% Test for Multivariate ARCH
% Input:
% res = TxK residual matrix
% lags = number of lagged ARCH effects
% K = number of variables in the underlying VAR/VEC model
%
% Function produces a table with results
%
% Reference: Lutkepohl (2005), New Introduction to Multiple Time Series
% Analysis, Springer. pp. 576-577
% Doornik and Hendry (1997), Modelling Dynamic Systems Using PcFiml 9.0 for
% Windows, International Thomson Business Press.
%
% Michael Bergman, 2018
%
uhat = res-ones(length(res),1)*mean(res);
i = 1;
% Initialize matrices
UUT = zeros(K*(K+1)/2,1);
i=1;
while i <=length(uhat)
utut = uhat(i,:)'*uhat(i,:);
  % the vech operator
  tmp  = [];
  for ii=1:K
     tmp = [tmp; utut(ii:end,ii)];
  end
  UUT = [UUT tmp];
  i = i+1;
end 
UUT=UUT(:,2:length(UUT));
% Create matrices of regressors
Y = UUT(:,1+lags:length(UUT));
Z = ones(1,length(uhat)-lags); 
i = 1;
while i <= lags
  Z = [Z; UUT(:,1+lags-i:length(UUT)-i)];
  i = i+1;
end
A=Y*Z'*inv(Z*Z');
omega = (Y-A*Z)*(Y-A*Z)';
omega = omega/length(Y);
% Compute omega0
a=mean(Y');
omega0 = (Y-a')*(Y-a')';
omega0 = omega0/length(Y);
R2 = 1-2/(K*(K+1))*sum(diag(omega*inv(omega0)));
VARCHLM = length(Y)*K*(K+1)*R2/2 ;
dof = lags*K^2*(K+1)^2/4;
display('Tests for Multivariate ARCH');
test=[VARCHLM; 1-chi2cdf(VARCHLM,dof); dof];
disp(table(categorical({'test statistic:' ; 'p-value' ; 'degrees of freedom'}),test,'VariableNames',{'Test' 'Doornik_Hendry'}));
end