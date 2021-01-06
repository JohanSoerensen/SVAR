function q = restrictions(B0inv,res_cov, Theta1)
% Long run multiplier
Upsilon = Theta1 * B0inv; 

% Restrictions - custommise these
q = [ vec(B0inv*B0inv'-res_cov);...
      B0inv(2,2);...                   % Short-run restrictions (insert indices in () )
      Upsilon(:,2);...                 % Long-run restrictions (insert indices in () )
      Upsilon(:,3);...                 % Long-run restrictions (insert indices in () )
      ];
end

