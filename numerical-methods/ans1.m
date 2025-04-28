clc; clearvars; close all;

% Constants and initial values
A = 3; 
B = -2; 
C = 3;

x_0 = 5; 
y_0 = 7;
max_iterations = 100;  
tolerance = 1e-5;      

iter_vector = [x_0; y_0]; 
iteration_history = iter_vector;
f_history = [];
          
syms x y;
f = A*x^2 - B*x*y + C*y^2 + x - y; 

% Create surface plot of the function
figure;
[X, Y] = meshgrid(-2:0.1:8, -2:0.1:8);
Z = A*X.^2 - B*X.*Y + C*Y.^2 + X - Y;
surf(X, Y, Z, 'DisplayName', 'Function Surface');
hold on;

% Calculate gradient
grad_f = gradient(f, [x, y]);

% Initialize table data
iterations = 0:max_iterations;
table_data = cell(max_iterations+1, 7);  % Changed to 7 columns
table_data{1,1} = 0;
table_data{1,2} = sprintf('(%.5f, %.5f)', x_0, y_0);
initial_grad = double(subs(grad_f, [x, y], [x_0; y_0]'));
table_data{1,3} = sprintf('(%.5f, %.5f)', initial_grad(1), initial_grad(2));
initial_grad_norm = initial_grad / norm(initial_grad);
table_data{1,4} = sprintf('(%.5f, %.5f)', initial_grad_norm(1), initial_grad_norm(2));
table_data{1,5} = '-';
table_data{1,6} = '-';
table_data{1,7} = '-';

% Main iteration loop
for i = 1:max_iterations
    % Calculate gradient and normalize
    grad_f = gradient(f, [x, y]);
    grad_val = subs(grad_f, [x, y], iter_vector');
    grad_norm = grad_val / norm(grad_val);
    
    % Calculate optimal alpha symbolically
    syms alpha
    x_next = iter_vector - alpha * grad_norm;
    
    % Substitute x_next into original function to get f(alpha)
    f_alpha = subs(f, [x, y], x_next');
    
    % Find optimal alpha
    df_dalpha = diff(f_alpha, alpha);
    alpha_optimal = solve(df_dalpha == 0, alpha);
    
    % Calculate next point and function value
    x_next = double(iter_vector - alpha_optimal * grad_norm);
    f_next = double(subs(f, [x, y], x_next'));
    
    % Store table data
    table_data{i+1,1} = i;
    table_data{i+1,2} = sprintf('(%.5f, %.5f)', iter_vector(1), iter_vector(2));
    grad_val_double = double(grad_val);
    table_data{i+1,3} = sprintf('(%.5f, %.5f)', grad_val_double(1), grad_val_double(2));
    grad_norm_double = double(grad_norm);
    table_data{i+1,4} = sprintf('(%.5f, %.5f)', grad_norm_double(1), grad_norm_double(2));
    table_data{i+1,5} = sprintf('%.5f', double(alpha_optimal));
    table_data{i+1,6} = sprintf('(%.5f, %.5f)', x_next(1), x_next(2));
    table_data{i+1,7} = sprintf('%.4e', norm(x_next - iter_vector));
    
    % Store history
    iteration_history = [iteration_history, x_next];
    f_history = [f_history, f_next];
    
    % Check convergence
    if i > 1
        position_change = norm(x_next - iter_vector);
        if position_change < tolerance
            fprintf('Converged after %d iterations\n', i);
            break;
        end
    end
    
    % Update for next iteration
    iter_vector = x_next;
end

% Display the table
fprintf('\nOptimization Progress:\n');
headers = {'Iteration', '(x,y)', '∇f(x,y)', '∇f(x,y) normalized', 'α', 'next(x,y)', 'Error'};
table_data = table_data(1:i+1,:);  % Trim unused rows
disp(array2table(table_data, 'VariableNames', headers));

% Plot initial point
z_initial = double(subs(f, [x, y], iteration_history(:,1)'));
plot3(iteration_history(1,1), iteration_history(2,1), z_initial, 'ro', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Initial Point');

% Plot final point
z_final = double(subs(f, [x, y], iteration_history(:,end)'));
plot3(iteration_history(1,end), iteration_history(2,end), z_final, 'go', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Final Point');

% Customize the surface plot
title('Function Surface with Optimization Path');
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
colorbar;
legend('Location', 'best');
grid on;
view(-30, 30);  

% Contour plot with optimization path
figure;
contour(X, Y, Z, 50); hold on;          % Plot contour lines

% Plot optimization path
plot(iteration_history(1,:), iteration_history(2,:), 'b-', 'LineWidth', 2);

% Plot initial and final points
plot(iteration_history(1,1), iteration_history(2,1), 'ro', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Initial Point');
plot(iteration_history(1,end), iteration_history(2,end), 'go', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Final Point');

% Customize plot
xlabel('x'); ylabel('y');
title('Contour + Optimization Path');
grid on;
legend('Location', 'northwest');

% Display results with 5 decimal places
fprintf('\nSolution found:\n');
fprintf('x = %.5f\n', iter_vector(1));
fprintf('y = %.5f\n', iter_vector(2));
fprintf('f(x,y) = %.5f\n', f_history(end));

% Calculate gradient at solution to verify minimum
grad_at_solution = double(subs(grad_f, [x, y], iter_vector'));
fprintf('\nGradient magnitude at solution: %.5e\n', norm(grad_at_solution));
