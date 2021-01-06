function [IRF_bs,FEVD_bs] = bootstrap(nrep,h,K,T,P,y,A,res,res_cov,cons,trend,restriction)
%BOOTSTRAP Summary of this function goes here
% nrep = number of reps
% alpha = significance level
% h = forecast horizon
% K = number of variables
% P = number of lags
% restriction = 1 for cholesky, 2 for custom restrictions
% cons = 1 if constant included, 0 if not
% trend = 1 if trend included, 0 if not

% preallocate
IRF_bs = nan(h+1,K*K,nrep);
FEVD_bs= nan(h,K*K,nrep);

%save constant settings for later modelling
cons_bs = cons;
trend_bs = trend;

%Calculate constants used in boostrapping sample construction
if cons == 1
cons = A(K*P+1,:)';  
end
if trend == 1
    tr = A(K*P+2,:)';
else 
    tr = trend;
end

% initialise
nvalid = 0;
ntrials = 0;

while nvalid<nrep
    % iterate total trials
    ntrials = ntrials+1;
    
    % print iteration number if restriction == 2
    if restriction == 2
        disp('This is trial number');
        disp(ntrials);
    end
    
    % break if no valid trials
    if nvalid==0 && ntrials >5000
        error('Max number of trials reached');
    end
    
    % draw boostrap residuals from estimated residuals
    idx = randi([1,T-P],T-P,1);                                             % index of random integers from [1;T] interval
    res_bs = res(idx,:);                                                    % reshuffle the estimated residuals 
    res_bs = [zeros(P,K);res_bs];                                           % add residuals (assumed zero) for inital conditions
    
    % de-mean boostrap residuals if no cons
    if cons == zeros(K,1)
        res_bs = res_bs - mean(res_bs);
    end 
    
    % simulate boostrap sample
    y_bs = NaN(T,K);                                                         % preallocate memory using trend and cons
    y_bs(1:P,:) = y(1:P,:);                                                  % fill in initial values from original sample
    
    % create bootstrap sample smart way
%     for t = P+1:T
%         for i=1:P
%         y_bs(t,:) = y_bs(t,:)+ y_bs(t-i,:)*A(1+(K*(i-1)):i*K,:);  % autoregressive terms
%         end
%         y_bs(t,:) = y_bs(t,:) + res_bs(t,:)+ cons' + t*tr'; % bs residuals + deterministics
%     end

%example below is for 3 lags
    for t = P+1:T
        % create bootstrap sample                                            
        y_bs(t,:) =  y_bs(t-1,:)*A(1:K,:) + y_bs(t-2,:)*A(K+1:2*K,:) + ...  % autoregressive terms
                    + y_bs(t-3,:)*A(2*K+1:3*K,:) + res_bs(t,:)+ ...                                        % bs residuals
                    cons' + t*tr';                                          % deterministics
    end
    
    % estimate VAR on bootstrap sample

    [A_bs,~,~,~,~,res_cov_bs,~,~,~] = mols(y_bs,P,cons_bs,trend_bs);
    
    % companion matrix
    Abs_comp = comp(A_bs,P);
    
    % check stability
    if max(eig(Abs_comp))>=1
        continue
    end
    
    % identification matrix
    switch restriction
        case 1
            B0inv_bs = chol(res_cov_bs)';
         
        case 2
            % Compute theta1 (theta1*B0inv gives the long run multipliers)
            Asumbs = Abs_comp(1:K,1:K);
            for i = 1:P-1
                Asumbs = Asumbs + Abs_comp(1:K,i*K+1:i*K+K);
            end
            theta1 = inv(eye(K)-Asumbs);
            
            % Function handle with long and short run restrictions
            restrictonsbs = @(B0inv_bs) [vec(B0inv_bs*B0inv_bs'-res_cov_bs);...                        
                      B0inv_bs(1,2:5)'; B0inv_bs(2,3:5)'; B0inv_bs(3,4:5)';... % Short-run restrictions
                      theta1(4,4)*B0inv_bs(4,5)+theta1(4,5)*B0inv_bs(5,5)]';  % Long-run restrictions

            % compute identification matrix (options already defined)
            B0inv_bs = fsolve(restrictonsbs,randn(K,K), options);
    end
    
    % check if identification is valid
    if max(max(abs(B0inv_bs*B0inv_bs'-res_cov_bs))) > 1e-5
        continue
    end
    
    %  count valid trials
    nvalid = nvalid+1;
    
    % bootstrap IRFs and FEVDs
    IRF_bs_draw = irfstruc(Abs_comp,B0inv_bs,K,P,h);                       % all IRFs
    IRF_bs(:,:,nvalid) = IRF_bs_draw(:,:);                                 % Collect IRF

    FEVD_bs_draw = fevd_irf(IRF_bs_draw,K,h);                              % All FEVDs
    FEVD_bs(:,:,nvalid) = FEVD_bs_draw(:,:);                           % share of variance in all vars attributed to mp shock
    
    % termination statement 
    if nvalid==nrep
        disp('------------------------------------------');
        disp('Great success! Loop ran without problems.');
        disp('Amount of valid trials reached:'); disp(nvalid);
        disp('Number of trials used:'); disp(ntrials);
        disp('------------------------------------------');
    end
end

% count invalid trials
ninvalid = ntrials - nvalid;
end

