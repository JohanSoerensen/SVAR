function [beta,se2_beta,se_beta,t_beta,res,res_cov,mAIC,mSIC,mHQC] = mols(y,p,cons,tr)
% Multivariate least squares estimator to estimate VAR models
%   Inputs:
%       y = time series matrix
%       p = preferred number of lags (int)
%       cons = dummy to include (1) or exclude (0) constant
%       tr = cons = dummy to include (1) or exclude (0) linear trend
%   Returns:
%       beta = parameter matrix
%       se_beta = standard errors of beta
%       t_beta = t ratios of beta
%       res_cov = variance-covariance matrix of residuals <-- SMALL SAMPLE
%       CORRECTION???
%       mAIC = multivariate Akaike information criterium
%       mSIC = multivariate Schwartz (/Bayesian) information criterium
%       mHQC = multivariate Hannan-Quinn information criterium
%   References: 
%       This function is based on the LSestimators.m file provided by Michael
%       Bergman on the Absalon course room. See also Kilian & Lutkepohl
%       (2017), ch. 2.

% a. dimensions
[~,K]   = size(y);
Kp      = K*p;

% b. dependent variable
dep     = y(p+1:length(y),:);


% c. independent variables
indep   = lagmatrix(y,1:p); 
indep   = indep(p+1:length(indep),:);

% d. deterministic terms
if cons == 1
    indep = [indep ones(length(indep),1)];                                  
end

if tr   == 1                                                                
    indep = [indep (1:length(dep))'];
end 

% Note: we need to find the number of parameters in each of the K equations
% in order to make the small sample correction of the covariance matrix
% we also need the number of obs
[T,Kp2] = size(indep);

% beta = inv(indep'*indep)*indep'*dep; which is equivalent to
beta    = indep\dep;

res     = dep-indep*beta;
res_cov = (res'*res)/(T-Kp2);
se2_beta = kron(res_cov,inv(indep'*indep));
se_beta = sqrt(diag(se2_beta));

if cons == 1 & tr == 0
t_beta  = reshape(beta, [(Kp+1)*K,1])./sqrt(diag(se2_beta));                 % vectorize beta to make consistent with covariance
end

if cons == 1 & tr == 1
t_beta  = reshape(beta, [(Kp+2)*K,1])./sqrt(diag(se2_beta));                 % vectorize beta to make consistent with covariance
end
    

% matrix in order to compute t-ratios
% multivariate information criteria - NB: assumes cons==1 & tr==0
logdet = log(det(res_cov)*(T-Kp)/(T-p));                                                 % auxillary variables
m = p;
mK2K   = m*K^2+K;

mAIC   = logdet + 2/T*mK2K;
mSIC   = logdet + log(T)/T*mK2K;
mHQC   = logdet + 2/T*(log(log(T)))*mK2K;

% disp('Standard LS')                                                         % print results
% display(beta','LS estimate')

end

