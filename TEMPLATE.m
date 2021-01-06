%% SVAR library file
%% 1 - Set up dataset
clc; clear all; 
clear; 
%Import y as matrix
% CHANGE THIS TO THE RIGHT y FILE
load ggdata.txt;                                                           
% Call the date in the file y
y = ggdata;
clear ggdata;
% remove rows with missing variables
y = rmmissing(y);

%Labels of variables and T, K dimensions
lbl = {'short term interest rates', 'medium term interest rates', ...
    'long term interest rates'};
% save dimension for later use
[T,K] = size(y);                                                         
%Define date vector
dt = datetime(1969,1,1,'format','QQQ-yyyy'):calmonths(1):datetime ...
    (1988,12,1,'format','QQQ-yyyy');
dt = dt';

%% 2 - Plot series (customize code depending on data)
%K different time series plotted: 
for i = 1:K
     % plot all variables by looping over
    subplot(K,1,i); 
    % matrix rows and titles 
    plot(dt,y(:,i),'Color','Black');                                        
    title(lbl(1,i))
    ylabel('pct.');
end

%Plot all series in one:
figure;
plot(dt,y);
legend(lbl);
title('Interest rates');
ylabel('pct.');

%Summary stats: use for loop to fill summary matrix
for i = 1:K                                                                
    summary(1,i) = mean(y(:,i));
    summary(2,i) = median(y(:,i));
    summary(3,i) = min(y(:,i));
    summary(4,i) = max(y(:,i));
    summary(5,i) = sum(y(:,i));
    summary(6,i) = std(y(:,i));
end 

stats = ["mean","median","min","max","obs","std"];                          % stringarray of relevant statistics
sum_stats = array2table(summary,'VariableNames',lbl);                       % create table from matrix 
sum_stats.Properties.RowNames = stats;                                      % statistics and variables
disp(sum_stats);

clear stats sum_stats summary

%% 3 - Unit root tests

%% 4 - Model specification - choosing number of lags
% Reasonable max number of lag lengths
pmax = 24;                                                                   
% const and or trend
cons = 1; 
tr = 0; 

% containers for information criteria
mAICs = NaN(pmax,1);
mSICs = NaN(pmax,1);
mHQCs = NaN(pmax,1);

% run models with various lag lengths
for p=1:pmax
    % Specify y length to make information criteria comparable across
    % samples, as all models run on same sample. 
    yy = y(pmax+1:end,:); 
	yy_lag = lagmatrix(y,1:p);
    yy_lag = yy_lag(pmax+1:end,:);
    
    % run multivariate OLS with constant and no trend. 
    % y / yy_lag / lag length / constant / trend
    [~,~,~,~,~,~,mAIC,mSIC,mHQC] = mols_lag_selection(yy,yy_lag,p,cons,tr);
    
    % append information criteria
    mAICs(p) = mAIC;
    mSICs(p) = mSIC;
    mHQCs(p) = mHQC;
end 
clear yy yy_lag

% table with information criteria
tab = array2table([(1:pmax)',mAICs, mSICs, mHQCs]);
tab.Properties.VariableNames = {'Lag', 'AIC', 'SIC', 'HQC'};
disp(tab);

% find preferred lag length
[~, pref_mAIC] = min(mAICs);
[~, pref_mSIC] = min(mSICs);
[~, pref_mHQC] = min(mHQCs);

% table of preferred lag lengths
tab = array2table([pref_mAIC, pref_mSIC, pref_mHQC]);
tab.Properties.VariableNames = {'AIC', 'BIC', 'HQC'};
tab.Properties.RowNames = ["Pref. lag"];
disp(tab);

%% 5 - lag length choice 
%Choose based on discretion and above tests
P = 3;

% estimate model with preferred lag length
[beta,se2_beta,se_beta,t_beta,res,res_cov,mAIC,mSIC,mHQC] = mols(y,P,cons,tr);

% reshape std. errors and T-ratios
se = reshape(se_beta,K*P+cons+tr,K);
se = se';
se = reshape(se,K*(K*P+cons+tr),1);

t = reshape(t_beta,K*P+cons+tr,K);
t = t';
t = reshape(t,K*(K*P+cons+tr),1);

%Print output 
output = array2table([vec(beta') se t]);
colnames = ["ParameterEstimates","StdErrors","TValues"];
output.Properties.VariableNames = colnames;
disp(output);

%plot and exclude P first obs from dt to make length equivalent
plot(dt((1+P):end),res); 
legend(lbl);

% consistency check with matlab built-in
Mdl = varm(K,P);
%outcomment if no trend 
%Mdl.Trend = NaN; 
EstMdl = estimate(Mdl,[y]);
summarize(EstMdl)
clear Mdl EstMdl Mdl.Trend 

%% 6 - Test for autocorrelation
portmanteus = NaN(K,2);                                                     % preallocate
max = 36;

% run for different lags
for i=P:max
    [stat,pval]=portmanteu(res,i,P);
    portmanteus(i,1) = stat;
    portmanteus(i,2) = pval;
end

% --- we need p values above 0.05!
tab = array2table([(1:max)' portmanteus]);                                 % create table with output
tab.Properties.VariableNames = {'j', 'TestStatistic', 'PValue'};
disp(tab);

%% 7 - Test for normality and heteroskedasticity
% Test for normality --- we need p values above 0.05!
multnorm(res);
% Test for heteroskedasticity --- we need p values above 0.05!
march(res,P,K);

%% 8 - Plotting of residuals to understand potential misspecification
%rows in figure
rows = 2;
%columns in figure
columns = K; 

figure; 
for i = 1:K
    hold on
    subplot(rows,columns,i)
    histfit(res(:,i))
    title(lbl(i))
end
for i = 1:K
    hold on
    subplot(rows,columns,i+K)
    qqplot(res(:,i))
end
clear columns rows

%% 9 - Test for cointegration
% intercept and linear trend in cointegration vector, linear trend in data (in levels)
[~,pValue,stat,cValue,mles] = jcitest(y,'model','H*','lags',P-1); 

% r = rank inferred by the jcitest above - 
% use economic theory if working with rates etc. Used in the Johansen
% constraint test before choosing final VEC model and final rank.
r = 2; 

%% 10 - Johansen constraint tests to test for linear intercept and trend
% Johansen constraint test - Bcon (test for removing linear trend such that d0=0)
% [h,pValue,stat,cValue,mles] = jcontest(Y,r,test,Cons,Name,Value)
% Tests:
% 'ACon'	Test linear constraints on A.
% 'AVec'	Test specific vectors in A.
% 'BCon'	Test linear constraints on B.
% 'BVec'	Test specific vectors in B.
[~,pValue,stat,cValue,mles] = jcontest(y,r,'Bcon',{[zeros(1,K) 1]'},...
                                            'lags',P-1,'model','H*');
display(mles.paramVals.B,'Restricted coint vec');
display(mles.paramVals.d0,'Estimated trend under restriction');
display([stat pValue],'LR-test & Pvalue (if above 0.05, no trend)');

% intercept in cointegration vector, linear trend in data (in levels)
[~,pValue,stat,cValue,mles] = jcitest(y,'model','H1','lags',P-1); 
% >>> Change rX depending on choice of rank <<< 
rLL_H1 = mles.r2.rLL;

% Intercept in cointegration vector --> variables stationary around a constant mean
% with no linear trend in data (in levels)
[~,pValue,stat,cValue,mles] = jcitest(y,'model','H1*','lags',P-1);
% >>> Change rX depending on choice of rank <<< 
rLL_H1star = mles.r2.rLL;

% LR test for c1 = 0 (no linear trend in data in levels if P > 0.05)
lrtest = [-2*(rLL_H1star - rLL_H1) chi2cdf(rLL_H1star - rLL_H1,K,'upper')];
lrtab = array2table(lrtest,'VariableNames',{'TestStat' 'PValue - if above 0.05 no trend in levels'});
display(lrtab);

% plot cointegration vector of preferred model (NB: adjust mles.rX)
B = mles.r2.paramVals.B;                                                  
c0 = mles.r2.paramVals.c0;

lbls = ["cointegration relation 1","cointegration relation 2"];             % create labels
plot(dt,y*B+repmat(c0',T,1));                                               % plot vectors
legend(lbls,'Location','NorthWest');                                        % add labels as legends
grid on   

% Choose preferred number of cointegrating vectors
r = 2;
% Choose preferred model
model = 'H1*'; %H1* - There are intercepts in the cointegrated series and 
              %there are no deterministic trends in the levels of the data.

%% 11 - Exclusion and stationarity tests
% create cell array of restrictions depending on choice of model above
if strcmp(model,'H2') | strcmp(model,'H1') | strcmp(model,'H') 
    R = cell(1,K);
    I = eye(K);

    for i=1:K
        R{1,i} = I(:,i);
    end
else
    R = cell(1,K+1);
    I = eye(K+1);

    for i=1:K+1
        R{1,i} = I(:,i);
    end
end

% stationarity
[h,pValue,stat] = jcontest(y,r,'BVec',R,'model',model,'lags',P-1);
disp('Tests for stationarity');
display(h,'logicals');
display(pValue, 'p-values'); 
display(stat, 'Test statistics');

% exclusion
[h,pValue,stat,~,mles3] = jcontest(y,r,'BCon',R,'model',model,'lags',P-1);

disp('Tests for exclusion');
%display(mles(3).paramVals.B, 'restricted beta vectors');
display(h,'logicals');
display(pValue, 'p-values');
display(stat, 'Test statistics');

%% 12 - weak exogeneity test
% Recreate R matrix for K regardless of model, because the constant /
% linear trend obviously doesnt adjust (we test A matrix not B, see section
% 9 for test specifications)
R = cell(1,K);
I = eye(K);

for i=1:K
    R{1,i} = I(:,i);
end

[h,pValue,stat] = jcontest(y,r,'ACon',R,'model',model,'lags',P-1);

disp('Tests for weak exogeneity');
display(h,'logicals');
display(pValue, 'p-values'); 
display(stat, 'Test statistics');

%% Structural VAR modelling! Woo 
%Start by implementing RNG and numerical optimiser settings
rng(1);
warning off 
options = optimset('TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',1e+10,'MaxIter',2e5);

%% 13 - SVAR: Cholesky decompositon
% Macroeconomic intuition will be used from this point. 
% We first show how to use a cholesky decomposition to form a structural 
% VAR model, in oder to do impulse response functions. 

% In this example, we have a dataset of the short, medium and long run
% interest rates, in that order. We use a cholesky decomposition with the 
% chol function, which implies that the short run rate is independent of the 
% medium and long run rate, the medium rate is independent of the long run 
% rate, and the long run rate is a function of the shorter rates. Later, we
% argue for different long run restrictions. 

% Get companion matrix
[A] = comp(beta,P);
% Get B0inv for cholesky decomposition
B0inv = chol(res_cov)';

%% 14 - SVAR: Common trends model
% First, run the preferred model and get mles output
[~,~,~,~,mles] = jcitest(y,'model',model,'lags',P-1);

%extract parameters for preferred rX (choice of rank is 2 in this case)
beta = mles.r2.paramVals.B;
alpha = mles.r2.paramVals.A;
gamma1 = mles.r2.paramVals.B1;
gamma2 = mles.r2.paramVals.B2;
sigma = mles.r2.EstCov;

%Compute orthogonal complements
beta_perp = null(beta');
alpha_perp = null(alpha');

%Compute Xi 
gamma_sum = gamma1 + gamma2;
Xi = beta_perp*inv(alpha_perp'*(eye(K)-gamma_sum)*beta_perp)*alpha_perp';

% Compute B0inv. Use the restrictions function to impose r = 2 long-run 
% and (K − r) = 1 short-run restrictions.
B0inv = fsolve( @(B0inv) restrictions(B0inv, sigma, Xi) ,randn(K,K),options);

display(B0inv, 'B0inv identified with long- and short-run restriction');
% Check that B0 equals qxq zero matrix for robustness
display(B0inv*B0inv'-sigma,'Check that B0inv is correct, should be a zero matrix');
% Check zero restrictions on Epsilon
display(Xi*B0inv,'Check that imposed Upsilon restrictions are valid');
% compute the variance-covariance matrix of the identified structural shocks
display(B0inv^-1 * sigma * (B0inv^-1)',...
'Check that variance-covariance matrix equals the identy matrix');

%% 15 - SVAR: Long run restrictions
% As above, macroeconomic intuition is used. Be careful with long run
% restriction assumptions, IRF and FEVD estimates can be very biased.

% First, compute Theta1 where Theta1*B0inv gives the long run multipliers
Asum = A(1:K,1:K);
for i = 1:P-1
    Asum = Asum + A(1:K,i*K+1:i*K+K);
end
theta1 = inv(eye(K)-Asum);

% Next, customise the long run restrictions

%% 16 - SVAR: Boostrapping IRF and FEVD to get confidence intervals with the delta method
% Most of the work here is done in the boostrap function. 
% Open and adjust if needed (crimerider style)

% nrep = number of reps
nrep = 500;
% h = forecast horizon
h = 24;
% K = number of variables - should be predefined by now
% P = number of lags - should be predefined by now
% cons = 1 if constant included, 0 if not
% trend = 1 if trend included, 0 if not
trend = tr; 
% Restriction = 1 for cholesky, 2 for custom restrictions
restriction = 1;

[A,~,~,~,res,res_cov,~,~,~] = mols(y,P,cons,trend);
%VAR Cholesky decomposition example
% Get companion matrix
B0inv = chol(res_cov)';
A_comp = comp(A,P);

[IRF_bs,FEVD_bs] = bootstrap(nrep,h,K,T,P,y,A,res,res_cov,cons,trend,restriction);

IRF = irfstruc(A_comp,B0inv,K,P,h);                  % all IRFs
FEVD = fevd_irf(IRF,K,h+1);                              % All FEVDs

%Explanation of FEVD output:
%Row 1,2,3... time periods in forecast horizon
%Column 1 = share of FEV of variable 1 attributed to a shock to variable 1
%Column 2 = share of FEV of variable 2 attributed to variable 1
%Column K+1 = Share of FEV of variable 1 attributed to shock in variable 2 
%(should be 0 for a cholesky decompositioned SVAR)
%Column K+2 = Share of FEV of variable 2 attributed to shock to variable 2
%Column 3 = 3rd dimension containing each bootstrap draw

%% 17 - Computing and plotting IRFs with confidence bands (Delta method)
% set significance level
alpha = 0.05; 
crit = -norminv(alpha/2);

% compute confidence bands for IRFs
IRF_CIhi = nan(h+1,K^2); % container for upper CI
IRF_CIlo = nan(h+1,K^2); % container for lower CI

for j=1:K^2
     for i=1:h+1
         IRF_CIlo(i,j,:) = IRF(i,j) - crit*(nanstd(IRF_bs(i,j,:)));                                   
         IRF_CIhi(i,j,:) = IRF(i,j) + crit*(nanstd(IRF_bs(i,j,:)));
     end
end

% plot IRFs with confidence bands
% Row 1 = shocks from variable 1
% Column 1 = shocks to variable 1   
figure;
title_IRF = (1:K^2); % replace with relevant string vector
for i=1:K^2
    subplot(K,K,i);
    hold on
    title(title_IRF(i),'Interpreter','latex')
    plot(0:h,IRF(:,i),'k-','Linewidth',2);
    plot(0:h,IRF_CIlo(:,i),'k:');
    plot(0:h,IRF_CIhi(:,i),'k:');
end

%% 18 - Computing FEVD with confidence bands (Delta & Efron method) 
%alpha = 0.32;

FEVD_CIhi = nan(h+1,K^2); % container for upper CI
FEVD_CIlo = nan(h+1,K^2); % container for lower CI

% Delta method
for j=1:K^2
     for i=1:h
         FEVD_CIlo(i,j,:) = FEVD(i,j) - crit*(nanstd(FEVD_bs(i,j,:)));                                   
         FEVD_CIhi(i,j,:) = FEVD(i,j) + crit*(nanstd(FEVD_bs(i,j,:)));
     end
end

figure;
title_FEVD = (1:K^2); % replace with relevant string vector
for i=1:K^2
    subplot(K,K,i);
    hold on
    title(title_FEVD(i),'Interpreter','latex')
    plot(0:h,FEVD(:,i),'k-','Linewidth',2);
    plot(0:h,FEVD_CIlo(:,i),'k:');
    plot(0:h,FEVD_CIhi(:,i),'k:');
end

% Efron method 
FEVD_CI_Efron = 100*quantile(FEVDbs_MPshock,[alpha/2,1-alpha/2],3);

