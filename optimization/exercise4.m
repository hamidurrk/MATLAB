%% MATLAB 1

clear; clc;

f = @(x) x^5 + 8*x^4 + 4*x^3 - 2*x - 10;
f_prime = @(x) 5*x^4 + 32*x^3 + 12*x^2 - 2;      % First derivative
f_double_prime = @(x) 20*x^3 + 96*x^2 + 24*x;     % Second derivative

x0 = 0;              % Initial guess
tol = 1e-6;          % Tolerance
max_iter = 100;      % Maximum iterations


x = x0;
for i = 1:max_iter
    f_prime_val = f_prime(x);
    f_double_prime_val = f_double_prime(x);
    
    if abs(f_double_prime_val) < 1e-12
        fprintf('Second derivative too small. Method failed.\n');
        break;
    end
    
    x_new = x - f_prime_val / f_double_prime_val;
    
    fprintf('Iteration %2d: x = %.10f, f''(x) = %.10e\n', i, x_new, f_prime(x_new));
    
    if abs(x_new - x) < tol
        fprintf('\nConverged!\n');
        break;
    end
    
    x = x_new;
end

fprintf('Critical point: x = %.10f\n', x);
fprintf('f(x) = %.10f\n', f(x));
fprintf('f''(x) = %.10e\n', f_prime(x));

figure;
subplot(2,1,1);
x_plot = linspace(-10, 2, 1000);
y_plot = arrayfun(f, x_plot);
plot(x_plot, y_plot, 'b-', 'LineWidth', 1.5);
hold on;
plot(x, f(x), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
grid on;
xlabel('x');
ylabel('f(x)');
title('Function f(x) = x^5 + 8x^4 + 4x^3 - 2x - 10');
legend('f(x)', 'Critical point');

subplot(2,1,2);
y_prime_plot = arrayfun(f_prime, x_plot);
plot(x_plot, y_prime_plot, 'b-', 'LineWidth', 1.5);
hold on;
plot(x, f_prime(x), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
yline(0, 'k--', 'LineWidth', 1);
grid on;
xlabel('x');
ylabel('f''(x)');
title('First Derivative f''(x)');
legend('f''(x)', 'Critical point', 'y=0');

%% MATLAB 2
clc; clear; close all;

% a)
function [x, z_max] = simplex_algorithm(A, b, c)
    [m, n] = size(A);
    
    if size(b, 2) > 1
        b = b';
    end
    if size(c, 2) > 1
        c = c';
    end
    
    if any(b < 0)
        error('All elements of b must be non-negative');
    end
    tableau = zeros(m + 1, n + m + 1);
    tableau(1:m, 1:n) = A;                    % Constraint coefficients
    tableau(1:m, n+1:n+m) = eye(m);           % Slack variables
    tableau(1:m, end) = b;                    % RHS
    tableau(end, 1:n) = -c';                  % Objective coefficients (negative for max)
    
    basic_vars = n + (1:m);
    
    tol = 1e-10;
    
    max_iter = 1000;
    iter = 0;
    
    while iter < max_iter
        iter = iter + 1;
        
        obj_row = tableau(end, 1:end-1);
        if all(obj_row >= -tol)
            break; % Optimal solution found
        end
        
        [min_val, entering_col] = min(obj_row);
        if min_val >= -tol
            break; % Optimal
        end
        
        pivot_col = tableau(1:m, entering_col);
        ratios = inf(m, 1);
        
        for i = 1:m
            if pivot_col(i) > tol
                ratios(i) = tableau(i, end) / pivot_col(i);
            end
        end
        
        [min_ratio, leaving_row] = min(ratios);
        
        if isinf(min_ratio)
            error('Problem is unbounded');
        end
        
        pivot_element = tableau(leaving_row, entering_col);
        
        tableau(leaving_row, :) = tableau(leaving_row, :) / pivot_element;
        
        for i = 1:m+1
            if i ~= leaving_row
                multiplier = tableau(i, entering_col);
                tableau(i, :) = tableau(i, :) - multiplier * tableau(leaving_row, :);
            end
        end
        
        basic_vars(leaving_row) = entering_col;
    end
    
    if iter >= max_iter
        warning('Maximum iterations reached');
    end
    
    x = zeros(n, 1);
    for i = 1:m
        if basic_vars(i) <= n
            x(basic_vars(i)) = tableau(i, end);
        end
    end
    
    z_max = tableau(end, end);
    
    fprintf('\nSimplex Algorithm Results:\n');
    fprintf('Optimal solution x:\n');
    disp(x);
    fprintf('Maximum value z = %.6f\n', z_max);
    fprintf('Number of iterations: %d\n\n', iter);
end

clc;
A = [1 1 0;
     0 0 2;
     0 1 1];
b = [10; 10; 5];
c = [3; 1; 2];

fprintf('Problem: Maximize z(x) = c^T * x\n');
fprintf('Subject to: A*x <= b, x >= 0\n\n');
fprintf('A = \n');
disp(A);
fprintf('b = \n');
disp(b);
fprintf('c = \n');
disp(c);

% b)
fprintf('(b) USING SIMPLEX_ALGORITHM()\n');
[x_simplex, z_simplex] = simplex_algorithm(A, b, c);

% c)
fprintf('(c) USING MATLAB''S LINPROG()\n');

f = -c;  % Negate for maximization
lb = zeros(3, 1);  % x >= 0

options = optimoptions('linprog', 'Display', 'off');
[x_linprog, z_linprog_neg] = linprog(f, A, b, [], [], lb, [], options);
z_linprog = -z_linprog_neg;  % Negate back to get maximum

fprintf('\nLinprog Results:\n');
fprintf('Optimal solution x:\n');
disp(x_linprog);
fprintf('Maximum value z = %.6f\n\n', z_linprog);

% Compare results
fprintf('         COMPARISON OF RESULTS\n');

fprintf('Solution comparison:\n');
fprintf('   Method              x1        x2        x3     Max z(x)\n');
fprintf('   --------------------------------------------------------\n');
fprintf('   Simplex:         %8.4f  %8.4f  %8.4f  %10.4f\n', ...
        x_simplex(1), x_simplex(2), x_simplex(3), z_simplex);
fprintf('   Linprog:         %8.4f  %8.4f  %8.4f  %10.4f\n', ...
        x_linprog(1), x_linprog(2), x_linprog(3), z_linprog);

fprintf('\nDifferences:\n');
diff_x = abs(x_simplex - x_linprog);
diff_z = abs(z_simplex - z_linprog);
fprintf('   |x_simplex - x_linprog| = [%.2e, %.2e, %.2e]\n', ...
        diff_x(1), diff_x(2), diff_x(3));
fprintf('   |z_simplex - z_linprog| = %.2e\n\n', diff_z);

if max(diff_x) < 1e-6 && diff_z < 1e-6
    fprintf('   RESULTS ARE THE SAME (within numerical tolerance)\n\n');
else
    fprintf('   RESULTS DIFFER\n\n');
end

fprintf('Constraint verification for simplex solution:\n');
Ax = A * x_simplex;
fprintf('   A*x = [%.4f, %.4f, %.4f]''\n', Ax(1), Ax(2), Ax(3));
fprintf('   b   = [%.4f, %.4f, %.4f]''\n', b(1), b(2), b(3));
if all(Ax <= b + 1e-6)
    fprintf('   All constraints satisfied\n\n');
else
    fprintf('   Constraints violated\n\n');
end
