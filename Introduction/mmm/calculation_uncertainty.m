time = [0, 1, 2, 3, 4, 5, 6]; 
A_conc = [5, 4.093654, 3.3516, 2.744058, 2.246645, 1.839397, 1.505971];
ln_A = log(A_conc);

p = polyfit(time, ln_A, 1);
slope = p(1);
intercept = p(2);
k = -p(1); 
n = length(time);
y_fit = polyval(p, time); 
residuals = ln_A - y_fit;
SSE = sum(residuals.^2); 
sigma2 = SSE / (n - 2); 
X_mean = mean(time);
SE_slope = sqrt(sigma2 / sum((time - X_mean).^2)); 

uncertainty_k = SE_slope;
disp(['Uncertainty in k: ±', num2str(uncertainty_k)]);
