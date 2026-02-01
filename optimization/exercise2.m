%% MATLAB 1
clc; clear;

% Profit coefficients (negative for minimization)
f = [-70; -40; -20];

% Inequality constraints A*x <= b
A = [ 8    6    1;
      4    2   0.4;
      2   1.5   1 ];

b = [100; 48; 16];

% Lower bounds (x >= 0)
lb = [0; 0; 0];

[x_opt, fval, exitflag] = linprog(f, A, b, [], [], lb);

x_opt;
max_profit = -fval

%% MATLAB 2
clear; clc;

function y_mode = lognormal_mode(mu, sigma)
    % Log-normal pdf
    pdf_ln = @(y) (1./(y*sigma*sqrt(2*pi))) .* ...
                  exp(-(log(y)-mu).^2 ./ (2*sigma^2));
    
    obj = @(y) -pdf_ln(y);
    
    y0 = exp(mu);
    
    obj_safe = @(y) obj(abs(y));
    
    y_mode = fminsearch(obj_safe, y0);
    y_mode = abs(y_mode);
end

cases = [
    0  1;   % a
    3  1;   % b
    0  2    % c
];

figure;

for i = 1:size(cases,1)
    mu = cases(i,1);
    sigma = cases(i,2);

    y_mode = lognormal_mode(mu, sigma);

    y_exact = exp(mu - sigma^2);

    fprintf('Case %d: mu = %.1f, sigma = %.1f\n', i, mu, sigma);
    fprintf('  Numerical mode = %.6f\n', y_mode);
    fprintf('  Exact mode     = %.6f\n\n', y_exact);

    y = linspace(0.001, exp(mu+3*sigma), 1000);
    pdf = (1./(y*sigma*sqrt(2*pi))) .* ...
          exp(-(log(y)-mu).^2 ./ (2*sigma^2));

    subplot(1,3,i)
    plot(y, pdf, 'LineWidth', 2)
    hold on
    plot(y_mode, ...
         (1/(y_mode*sigma*sqrt(2*pi))) * ...
         exp(-(log(y_mode)-mu)^2/(2*sigma^2)), ...
         'ro', 'MarkerSize', 8, 'LineWidth', 2)
    hold off

    title(sprintf('\\mu = %.1f, \\sigma = %.1f', mu, sigma))
    xlabel('y')
    ylabel('pdf')
    grid on
end
