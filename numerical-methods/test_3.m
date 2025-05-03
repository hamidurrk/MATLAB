clc;

function [x_min, y_min, iterations] = newton_analytical(A, B, C, x0, y0, tol, max_iter)
    x = x0;
    y = y0;
    for iterations = 1:max_iter
        % Analytical gradient
        g_x = 2*A*x - B*y + 1;
        g_y = -B*x + 2*C*y - 1;
        
        % Optimal alpha calculation
        numerator = g_x^2 + g_y^2;
        denominator = 2*A*g_x^2 - 2*B*g_x*g_y + 2*C*g_y^2;
        alpha = numerator / denominator;
        
        % Update coordinates
        x_new = x - alpha * g_x;
        y_new = y - alpha * g_y;
        
        % Check convergence
        if abs(x_new - x) < tol && abs(y_new - y) < tol
            break;
        end
        x = x_new;
        y = y_new;
    end
    x_min = x;
    y_min = y;
end

function [x_min, y_min, iterations] = newton_numerical(f, A, B, C, x0, y0, tol, max_iter, h)
    x = x0;
    y = y0;
    for iterations = 1:max_iter
        % Numerical gradient (forward difference)
        f_current = f(A, B, C, x, y);
        g_x = (f(A, B, C, x + h, y) - f_current) / h;
        g_y = (f(A, B, C, x, y + h) - f_current) / h;

        % Function to find alpha where dot product is zero
        fun = @(alpha) compute_g_alpha(alpha, x, y, g_x, g_y, f, A, B, C, h);

        % Find alpha using fzero
        options = optimset('TolX', 1e-8);
        try
            alpha = fzero(fun, 0, options);
        catch
            alpha = 0.1; % Fallback if root not found
        end

        % Update coordinates
        x_new = x - alpha * g_x;
        y_new = y - alpha * g_y;

        % Check convergence
        if abs(x_new - x) < tol && abs(y_new - y) < tol
            break;
        end
        x = x_new;
        y = y_new;
    end
    x_min = x;
    y_min = y;
end

function g_alpha = compute_g_alpha(alpha, x, y, g_x, g_y, f, A, B, C, h)
    x_new = x - alpha * g_x;
    y_new = y - alpha * g_y;

    % Numerical gradient at new point
    f_new = f(A, B, C, x_new, y_new);
    g_x_new = (f(A, B, C, x_new + h, y_new) - f_new) / h;
    g_y_new = (f(A, B, C, x_new, y_new + h) - f_new) / h;

    % Dot product with old gradient
    g_alpha = g_x_new * g_x + g_y_new * g_y;
end

% Define the function f
function val = f(A, B, C, x, y)
    val = A*x^2 - B*x*y + C*y^2 + x - y;
end

% Constants and initial values
A = 3; B = -2; C = 3;
x0 = 5; y0 = 7;
tol = 1e-5;
max_iter = 1000;
h = 1e-20;         % Initial guess
x_target = -0.25;          % Known minimum x-coordinate
y_target = 0.25; 

% Part 1: Analytical
[x_analytical, y_analytical, iter_analytical] = newton_analytical(A, B, C, x0, y0, tol, max_iter);

% Part 2: Numerical
% [x_numerical, y_numerical, iter_numerical] = newton_numerical(@f, A, B, C, x0, y0, tol, max_iter, h);

disp(['Analytical Minimum: (', num2str(x_analytical), ', ', num2str(y_analytical), ') in ', num2str(iter_analytical), ' iterations']);
% disp(['Numerical Minimum:  (', num2str(x_numerical), ', ', num2str(y_numerical), ') in ', num2str(iter_numerical), ' iterations']);


% Test different h values
h_values = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14];
iterations = zeros(size(h_values));
x_min_values = zeros(size(h_values));
y_min_values = zeros(size(h_values));

for i = 1:length(h_values)
    h = h_values(i);
    [x_numerical, y_numerical, iter] = newton_numerical(@f, A, B, C, x0, y0, tol, max_iter, h);
    iterations(i) = iter;
    x_min_values(i) = x_numerical;
    y_min_values(i) = y_numerical;
    fprintf('h = %e → Iterations = %d\n', h, iter);
    disp(['Numerical Minimum:  (', num2str(x_numerical), ', ', num2str(y_numerical), ') in ', num2str(iter), ' iterations']);
end

% Compute absolute errors
x_errors = abs(x_min_values - x_target);
y_errors = abs(y_min_values - y_target);
figure;

% Subplot 1
subplot(2, 1, 1);
loglog(h_values, x_errors, '-o'); % Use loglog for both axes
xlabel('h (log scale)');
ylabel('|x_{numerical} - x_{target}| (log scale)');
title('Logarithmic Absolute Error in x_{numerical} vs h');
yline(1e-5, '--', 'Threshold', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom'); % Add dashed line

% Subplot 2
subplot(2, 1, 2);
loglog(h_values, y_errors, '-o'); % Use loglog for both axes
xlabel('h (log scale)');
ylabel('|y_{numerical} - y_{target}| (log scale)');
title('Logarithmic Absolute Error in y_{numerical} vs h');
yline(1e-5, '--', 'Threshold', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom'); % Add dashed line

% Plot Number of Iterations vs h (Discrete Plot)
figure;
semilogx(h_values, iterations, '-o', 'MarkerSize', 8, 'LineWidth', 1.5); % Use semilogx for discrete points
xlabel('h (log scale)');
ylabel('Number of Iterations');
title('Number of Iterations vs h');
grid on;

% Adjust y-axis range for better readability
ylim([0, max(iterations) + 4]); % Extend the y-axis range slightly above the maximum value

% Highlight the area where absolute error exceeds 10^-5
% hold on;
% threshold_indices = find(x_errors > 1e-5 | y_errors > 1e-5); % Find indices where error exceeds threshold
% if ~isempty(threshold_indices)
%     for idx = threshold_indices
%         xline(h_values(idx), '--r', 'LineWidth', 1.5); % Add vertical red dashed lines
%     end
% end
legend('Iterations', 'Error > 10^{-5}');
hold off;