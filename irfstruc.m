function [irf]=irfstruc(A,B0inv,K,p,h)
% Process to invert VAR to VMA
% Note: A is the companion matrix
%       K = number of var's
%       p = number of lags
%       h = horizon
if p==1;
  jmat = eye(K);
end;
if p>1;
jmat = [eye(K,K) zeros(K,K*(p-1))];
end;
irf = reshape( jmat*(A^0)*jmat'*B0inv, [1,K*K]   );
for i=1:h
   irf = [ irf; reshape( jmat*(A^i)*jmat'*B0inv, [1,K*K]   )];
end
end