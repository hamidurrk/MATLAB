% filepath: c:\Users\hamid\OneDrive\Documents\MATLAB\numerical-methods\ans2.m
clc; clearvars; close all;

% ---- Problem data ----
A = 3; B = -2; C = 3;
f  = @(x,y) A*x.^2 - B*x.*y + C*y.^2 + x - y;  % objective function

% Values of h to test
h_values = logspace(-1, -8, 8);  % 10 values from 1e-1 to 1e-10
final_xy = zeros(2, length(h_values));  % To store final x, y for each h

for h_idx = 1:length(h_values)
    h = h_values(h_idx);  % Current h value

    % forward-difference approximations
    dfdx_fd = @(x,y) (f(x+h,y) - f(x,y)) / h;
    dfdy_fd = @(x,y) (f(x,y+h) - f(x,y)) / h;

    % ---- Newton with dynamic alpha ----
    max_iter  = 100;
    tol       = 1e-5;
    x = 5; y = 7;  % starting point

    for k = 1:max_iter
        %--- 1) build numerical gradient and normalize
        gx = dfdx_fd(x,y);
        gy = dfdy_fd(x,y);
        g = [gx; gy];
        d = g / norm(g);

        %--- 2) form phi'(alpha)
        phi_prime = @(alpha) ...
          dfdx_fd( x - alpha*d(1), y - alpha*d(2) ) * d(1) + ...
          dfdy_fd( x - alpha*d(1), y - alpha*d(2) ) * d(2);

        %--- 3) bracket until sign change for phi'
        a_lo = 0;  f_lo = phi_prime(a_lo);
        a_hi = 1;  f_hi = phi_prime(a_hi);
        while f_lo * f_hi > 0
            a_hi = 2*a_hi;
            f_hi = phi_prime(a_hi);
            if a_hi > 1e6
                error('Could not bracket phi'' root.');
            end
        end

        %--- 4) solve phi'(alpha) = 0
        alpha_opt = fzero(phi_prime, [a_lo, a_hi]);

        %--- 5) update (x,y)
        x_new = x - alpha_opt * d(1);
        y_new = y - alpha_opt * d(2);

        % check convergence
        if norm([x_new - x; y_new - y]) < tol
            break
        end

        % prepare for next iteration
        x = x_new;
        y = y_new;
    end

    % Store final x, y for this h
    final_xy(:, h_idx) = [x; y];
end

% ---- Plot results ----
figure;

% Plot x vs h
semilogx(h_values, final_xy(1,:), '-o', 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'DisplayName', 'Final x');
hold on;

% Plot y vs h
semilogx(h_values, final_xy(2,:), '-s', 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'DisplayName', 'Final y');

% Add labels, title, and legend
xlabel('h (log scale)', 'FontSize', 12);
ylabel('Final x, y', 'FontSize', 12);
title('Effect of h on Final x and y', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);

% Improve grid and formatting
grid on;
set(gca, 'FontSize', 12);