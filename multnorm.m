function [norm]=multnorm(umat);
%{
**
**  Purpose:   Multivariate test on normality.
**
**
**  H0:        gaussian data generation process
**  Usage:     [norm] = multnorm(umat);
**
**  Input:     umat is a TxK matrix of residuals
**
**  Output:    NormPort is a 7x2 vector (Hendry and Lutkepohl standardization) 
**             |value|p-value|df|skew|pskew|kurt|pkurt|
**
**
**  Reference: Lutkepohl (1993), Introduction to Multiple Time Series
**             Analysis, 2ed, p. 150. (Cholesky variant); 
**             Doornik/Hansen (1994) (Jordan variant)
**             
** Based on Jmulti Gauss procedure
** Michael Bergman
%}
  [n, k] = size(umat);
  umat = umat - mean(umat);
  Su = 1/n*(umat'*umat);
  [Pmat,lambda] = eig(Su);
  i = 1;
  x = sqrt(diag(Pmat'*Pmat));
  [ rP cP ]=size(Pmat);
  while i <= cP;
    Pmat(:,i) = Pmat(:,i)./x(i);
    i = i + 1;
  end  
  Q = Pmat*(lambda.^0.5)*Pmat';
  v1= inv(Q) * (umat)';
  v2= inv(chol(Su)') * (umat)';
  b21 = (sum(v1'.^4)/n)';
  b11 = (sum(v1'.^3)/n)';
  b22 = (sum(v2'.^4)/n)';
  b12 = (sum(v2'.^3)/n)';
  l11 = n*b11'*b11/6; 
  pskew1= 1-chi2cdf( l11,k ); 
  l12 = n*b12'*b12/6; 
  pskew2 = 1-chi2cdf( l12,k ); 
  l21 = n*(b21 - 3)'*(b21 - 3)/24; 
  pkurt1 = 1-chi2cdf( l21,k ); 
  l22 = n*(b22 - 3)'*(b22 - 3)/24; 
  pkurt2 = 1-chi2cdf( l22,k ); 
  NormDf = 2*k;  
  l31 = l11 + l21;
  Normpv1 = 1-chi2cdf( l31,NormDf );
  l32 = l12 + l22;
  Normpv2 = 1-chi2cdf( l32,NormDf );
  norm = [ l31 l32; Normpv1 Normpv2; NormDf NormDf; l11 l12; pskew1 pskew2; l21 l22; pkurt1 pkurt2];
  display('Tests for non-normality');
  disp(table(categorical({'joint test statistic:' ; 'p-value' ; 'degrees of freedom' ; 'Skewness only' ; 'p-value' ; 'kurtosis only'; 'p-value'}),norm(:,1),norm(:,2),'VariableNames',{'Test' 'Doornik_Hansen' 'Lutkepohl'}));
end
