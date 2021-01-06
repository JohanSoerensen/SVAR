function [A]=comp(Beta,p);
% This function computes the companion matrix A using the estimated
% A1, A2, A3,...,Ap matrices.
% Note: It is not necessary to remove coefficients associated
% with constant and linear trend terms given that the input matrix
% Beta = [A1, A2, A3,...,Ap, Constant, Trend, Exogeneous ]'.
%
% Michael Bergman December 2018
%
[K,Kp]=size(Beta');  % loading Beta' so we need to take transpose
A = [ Beta(1:K*p,1:K)'; [ eye((p-1)*K) zeros((p-1)*K,K) ]];
end
