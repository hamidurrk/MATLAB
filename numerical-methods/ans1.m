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

% Display results with 5 decimal places
fprintf('\nSolution found:\n');
fprintf('x = %.5f\n', iter_vector(1));
fprintf('y = %.5f\n', iter_vector(2));
fprintf('f(x,y) = %.5f\n', f_history(end));

% Calculate gradient at solution to verify minimum
grad_at_solution = double(subs(grad_f, [x, y], iter_vector'));
fprintf('\nGradient magnitude at solution: %.5e\n', norm(grad_at_solution));

% Plot convergence history
figure;
semilogy(1:length(f_history), abs(diff([f_history(1), f_history])), 'b.-');
title('Convergence History (Log Scale)');
xlabel('Iteration');
ylabel('Change in Function Value');
grid on;
