%% Task 1
clear; clc; close all;

T = readtable('Advertising.csv');

% first column is an index column
if any(strcmpi(T.Properties.VariableNames{1}, {'Var1','X','Index'})) && width(T) > 4
    T(:,1) = [];
end

% extracting the variables
TV = T.TV;
Radio = T.radio;
Newspaper = T.newspaper;
Sales = T.sales;

n = length(Sales);
X = [ones(n,1), TV, Radio, Newspaper];
y = Sales;
p = size(X,2);   % number of parameters including intercept

disp('Data loaded successfully.');
disp(['Number of observations: ', num2str(n)]);
disp(['Number of parameters (including intercept): ', num2str(p)]);
disp(' ');

% TASK 1 starts here
disp('TASK 1');

% OLS estimates
beta_OLS = (X' * X) \ (X' * y);
yhat_OLS = X * beta_OLS;
residuals_OLS = y - yhat_OLS;

% sums of squares
SSE = sum(residuals_OLS.^2);                  % error 
SST = sum((y - mean(y)).^2);                  % total 
SSR = sum((yhat_OLS - mean(y)).^2);           % regression 

% variance estimate
sigma2_hat = SSE / (n - p);
sigma_hat = sqrt(sigma2_hat);

% cov matrix and standard errors
Cov_beta_OLS = sigma2_hat * inv(X' * X);
SE_beta_OLS = sqrt(diag(Cov_beta_OLS));

% t critical value for 95% CI
alpha = 0.05;
tcrit = tinv(1 - alpha/2, n - p);

% CIs
CI_OLS = [beta_OLS - tcrit * SE_beta_OLS, beta_OLS + tcrit * SE_beta_OLS];

R2 = 1 - SSE / SST;
adjR2 = 1 - (SSE/(n-p)) / (SST/(n-1));

disp('OLS coefficient estimates:');
paramNames = {'Intercept'; 'TV'; 'Radio'; 'Newspaper'};
for i = 1:p
    disp([paramNames{i}, ': ', num2str(beta_OLS(i), '%.6f')]);
end
disp(' ');

disp('OLS standard errors:');
for i = 1:p
    disp([paramNames{i}, ': ', num2str(SE_beta_OLS(i), '%.6f')]);
end
disp(' ');

disp('95% confidence intervals for OLS coefficients:');
for i = 1:p
    disp([paramNames{i}, ': [', num2str(CI_OLS(i,1), '%.6f'), ', ', ...
          num2str(CI_OLS(i,2), '%.6f'), ']']);
end
disp(' ');

% Hypothesis test: H0: beta1 = beta2 = beta3 = 0
% Restricted model: y = beta0 + error
X_restricted = ones(n,1);
beta_restricted = (X_restricted' * X_restricted) \ (X_restricted' * y);
yhat_restricted = X_restricted * beta_restricted;
SSE_restricted = sum((y - yhat_restricted).^2);

q = 3;  % number of restrictions
F_stat = ((SSE_restricted - SSE)/q) / (SSE/(n-p));
p_value_F = 1 - fcdf(F_stat, q, n-p);

disp('Joint hypothesis test: H0: beta_TV = beta_Radio = beta_Newspaper = 0');
disp(['F-statistic = ', num2str(F_stat, '%.6f')]);
disp(['p-value     = ', num2str(p_value_F, '%.6e')]);

if p_value_F < alpha
    disp('Decision: Reject H0 at alpha = 0.05');
else
    disp('Decision: Fail to reject H0 at alpha = 0.05');
end
disp(' ');

disp('Model fit statistics (OLS):');
disp(['SSE          = ', num2str(SSE, '%.6f')]);
disp(['sigma^2_hat  = ', num2str(sigma2_hat, '%.6f')]);
disp(['R^2          = ', num2str(R2, '%.6f')]);
disp(['Adjusted R^2 = ', num2str(adjR2, '%.6f')]);
disp(' ');

% QQ plot of residuals
figure;
qqplot(residuals_OLS);
title('QQ-Plot of OLS Residuals');
grid on;

%% TASK 2: BAYESIAN LINEAR REGRESSION
disp('TASK 2: BAYESIAN LINEAR REGRESSION');
%Using a conjugate Normal-Inverse-Gamma prior.
%This prior does NOT recover the OLS case because it is proper and finite.

% Model:
% y | beta, sigma^2 ~ N(X*beta, sigma^2 I)
% beta | sigma^2 ~ N(b0, sigma^2 * V0)
% sigma^2 ~ Inv-Gamma(a0, d0)

b0 = [0; 0; 0; 0];                    % prior mean
V0 = diag([100, 1, 1, 1]);            % finite prior covariance scaling
a0 = 2;                               % prior shape
d0 = 1;                               % prior scale

% Posterior updates
V0_inv = inv(V0);
Vn = inv(V0_inv + X' * X);
bn = Vn * (V0_inv * b0 + X' * y);
an = a0 + n/2;
dn = d0 + 0.5 * (y' * y + b0' * V0_inv * b0 - bn' * inv(Vn) * bn);

sigma2_post_mean = dn / (an - 1);

Cov_beta_post = (dn / (an - 1)) * Vn;

% Marginal posterior of beta is multivariate t
% beta_j | y follows Student-t with df = 2*an
df_post = 2 * an;
SE_beta_post = sqrt(diag((dn/an) * Vn));   % for t-based intervals

% 95% CIs
tcrit_post = tinv(1 - alpha/2, df_post);
CI_Bayes = [bn - tcrit_post * SE_beta_post, bn + tcrit_post * SE_beta_post];

yhat_Bayes = X * bn;
residuals_Bayes = y - yhat_Bayes;

disp('Posterior mean estimates of coefficients:');
for i = 1:p
    disp([paramNames{i}, ': ', num2str(bn(i), '%.6f')]);
end
disp(' ');

disp('Posterior standard deviations of coefficients:');
for i = 1:p
    disp([paramNames{i}, ': ', num2str(sqrt(Cov_beta_post(i,i)), '%.6f')]);
end
disp(' ');

disp('95% Bayesian credible intervals for coefficients:');
for i = 1:p
    disp([paramNames{i}, ': [', num2str(CI_Bayes(i,1), '%.6f'), ', ', ...
          num2str(CI_Bayes(i,2), '%.6f'), ']']);
end
disp(' ');

disp('Posterior summary for sigma^2:');
disp(['Posterior shape an = ', num2str(an, '%.6f')]);
disp(['Posterior scale dn = ', num2str(dn, '%.6f')]);
disp(['Posterior mean E[sigma^2 | y] = ', num2str(sigma2_post_mean, '%.6f')]);
disp(' ');

% Optional Bayesian fit measures
SSE_Bayes = sum((y - yhat_Bayes).^2);
R2_Bayes = 1 - SSE_Bayes / SST;
adjR2_Bayes = 1 - (SSE_Bayes/(n-p)) / (SST/(n-1));

disp('Fit based on posterior mean predictions:');
disp(['SSE          = ', num2str(SSE_Bayes, '%.6f')]);
disp(['R^2          = ', num2str(R2_Bayes, '%.6f')]);
disp(['Adjusted R^2 = ', num2str(adjR2_Bayes, '%.6f')]);
disp(' ');

% QQ plot of Bayesian residuals 
figure;
qqplot(residuals_Bayes);
title('QQ-Plot of Residuals from Bayesian Posterior Mean Fit');
grid on;

disp(' ');
disp(' ');
disp(' ');

disp('COMPARISON: OLS vs BAYESIAN POSTERIOR MEAN');

% Table 1: OLS results
OLS_Table = table(paramNames, beta_OLS, CI_OLS(:,1), CI_OLS(:,2), ...
    'VariableNames', {'Parameter','OLS_Estimate','CI_Lower','CI_Upper'});

disp('OLS results table:');
disp(OLS_Table);
disp(' ');

% Table 2: Bayesian results
Bayes_Table = table(paramNames, bn, CI_Bayes(:,1), CI_Bayes(:,2), ...
    'VariableNames', {'Parameter','Bayes_PostMean','CI_Lower','CI_Upper'});

disp('Bayesian results table:');
disp(Bayes_Table);
disp(' ');