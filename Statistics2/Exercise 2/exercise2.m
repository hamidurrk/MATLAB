%% Question 1
clear; clc; close all;

T = readtable('Exercise 2\basketball_data.csv');

players = unique(T.player);
m = length(players);

n = zeros(m,1);
ybar = zeros(m,1);
s2 = zeros(m,1);

% sample size, sample mean, and sample variance for each player
for i = 1:m
    yi = T.points(T.player == players(i));
    n(i) = length(yi);
    ybar(i) = mean(yi);
    s2(i) = var(yi);   
end

% methods-of-moments Empirical Bayes estimates

% pooled within-player variance estimate for sigma^2
sigma2_hat = sum((n - 1) .* s2) / sum(n - 1);

% mu from average of player means
mu_hat = mean(ybar);

% tau^2 from variance of player means
s2_ybar = var(ybar);   % sample variance across player means
tau2_hat = s2_ybar - sigma2_hat * mean(1 ./ n);

% preventing negative estimate
tau2_hat = max(tau2_hat, 0);

% posterior means
B = tau2_hat ./ (tau2_hat + sigma2_hat ./ n);
post_mean = B .* ybar + (1 - B) .* mu_hat;

Results = table(players, n, ybar, B, post_mean, ...
    'VariableNames', {'Player','Games','SampleMean','ShrinkageFactor','PosteriorMean'});

disp('Empirical Bayes estimates:')
fprintf('mu_hat     = %.4f\n', mu_hat);
fprintf('sigma2_hat = %.4f\n', sigma2_hat);
fprintf('tau2_hat   = %.4f\n\n', tau2_hat);

disp(Results);

%sample mean vs posterior mean
figure;
scatter(ybar, post_mean, 80, 'filled');
hold on;
plot([min(ybar)-1, max(ybar)+1], [min(ybar)-1, max(ybar)+1], 'r--', 'LineWidth', 1.5);
xlabel('Raw sample mean');
ylabel('Posterior mean (shrunk)');
title('Empirical Bayes shrinkage of player scoring means');
grid on;

% points with player ID
for i = 1:m
    text(ybar(i)+0.05, post_mean(i), sprintf('Player %d', players(i)));
end

% amount of shrinkage vs number of games
figure;
scatter(n, ybar - post_mean, 80, 'filled');
xlabel('Number of games');
ylabel('Raw mean - Posterior mean');
title('Shrinkage amount vs number of games');
grid on;

%% Question 2
clear; clc; close all;

alpha = 0.05;
cred_mass = 1 - alpha;

% Normal N(0.5, 2) ---------------------
mu = 0.5;
sigma = sqrt(2);

et_norm = norminv([alpha/2, 1-alpha/2], mu, sigma);

% For a symmetric unimodal normal, HPD = equal-tailed
hpd_norm = et_norm;

fprintf('Normal N(0.5, 2)\n');
fprintf('Equal-tailed CI: [%.4f, %.4f]\n', et_norm(1), et_norm(2));
fprintf('HPD CI:          [%.4f, %.4f]\n\n', hpd_norm(1), hpd_norm(2));

% Beta(5,2) ---------------------
aB = 5; bB = 2;
et_beta = betainv([alpha/2, 1-alpha/2], aB, bB);

% HPD by shortest 95% interval over quantiles
p = linspace(0, alpha, 200000);
lower_beta = betainv(p, aB, bB);
upper_beta = betainv(p + cred_mass, aB, bB);
[~, idxB] = min(upper_beta - lower_beta);
hpd_beta = [lower_beta(idxB), upper_beta(idxB)];

fprintf('Beta(5,2)\n');
fprintf('Equal-tailed CI: [%.4f, %.4f]\n', et_beta(1), et_beta(2));
fprintf('HPD CI:          [%.4f, %.4f]\n\n', hpd_beta(1), hpd_beta(2));

% Gamma(2,2) ---------------------
% shape = 2, scale = 2
aG = 2; bG = 2;
et_gamma = gaminv([alpha/2, 1-alpha/2], aG, bG);

lower_gamma = gaminv(p, aG, bG);
upper_gamma = gaminv(p + cred_mass, aG, bG);
[~, idxG] = min(upper_gamma - lower_gamma);
hpd_gamma = [lower_gamma(idxG), upper_gamma(idxG)];

fprintf('Gamma(2,2)  [shape=2, scale=2]\n');
fprintf('Equal-tailed CI: [%.4f, %.4f]\n', et_gamma(1), et_gamma(2));
fprintf('HPD CI:          [%.4f, %.4f]\n\n', hpd_gamma(1), hpd_gamma(2));

% plotting distributions and intervals 
x1 = linspace(mu - 4*sigma, mu + 4*sigma, 1000);
y1 = normpdf(x1, mu, sigma);

x2 = linspace(0, 1, 1000);
y2 = betapdf(x2, aB, bB);

x3 = linspace(0, gaminv(0.999, aG, bG), 1000);
y3 = gampdf(x3, aG, bG);

figure;
plot(x1, y1, 'LineWidth', 1.5); hold on;
xline(et_norm(1), '--r'); xline(et_norm(2), '--r');
xline(hpd_norm(1), ':k'); xline(hpd_norm(2), ':k');
title('Normal N(0.5,2)');
legend('pdf','Equal-tailed','Equal-tailed','HPD','HPD');
grid on;

figure;
plot(x2, y2, 'LineWidth', 1.5); hold on;
xline(et_beta(1), '--r'); xline(et_beta(2), '--r');
xline(hpd_beta(1), ':k'); xline(hpd_beta(2), ':k');
title('Beta(5,2)');
legend('pdf','Equal-tailed','Equal-tailed','HPD','HPD');
grid on;

figure;
plot(x3, y3, 'LineWidth', 1.5); hold on;
xline(et_gamma(1), '--r'); xline(et_gamma(2), '--r');
xline(hpd_gamma(1), ':k'); xline(hpd_gamma(2), ':k');
title('Gamma(2,2)');
legend('pdf','Equal-tailed','Equal-tailed','HPD','HPD');
grid on;

%% Comment on symmetry and overlap

% Normal distribution is symmetric, so the HPD and equal-tailed intervals are exactly the same.

% Beta(5,2) is left-skewed on [0,1], so the HPD and equal-tailed intervals are different, but still overlap a lot.

% Gamma(2,2) is strongly right-skewed, so the difference is more noticeable:

%      - equal-tailed interval stretches farther into the right tail

%      - HPD interval is shorter and stays more concentrated where density is highest