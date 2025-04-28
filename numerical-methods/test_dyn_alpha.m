clc; clearvars; 

% Constants and initial values
x_0 = 1; 
y_0 = 3;

iter_vector = [x_0; y_0]; 

% Function
syms x y; 
f = 25*x^2 + y^2;
grad_f = gradient(f, [x, y]); 
disp('Gradient of the function:');
disp(grad_f);

grad_val = subs(grad_f, [x, y], iter_vector'); 
grad_norm = grad_val / norm(grad_val);

disp('Gradient value at the initial point:');
disp(grad_val);
disp('Normalized gradient value:');
disp(double(grad_norm));

% syms alpha;
% x_next = double(iter_vector) - alpha * double(grad_norm);
% disp('Decimal representation of x_next:');
% disp(vpa(x_next, 10)); 
% iter_result = double(iter_vector - alpha * grad_norm)


% Calculate x(1) symbolically
syms alpha
grad_norm_double = double(grad_norm);
x_next = iter_vector - alpha * grad_norm_double;

% Substitute x_next into original function to get f(alpha)
f_alpha = subs(f, [x, y], x_next');

% Take derivative of f with respect to alpha and solve for alpha
df_dalpha = diff(f_alpha, alpha);
alpha_optimal = solve(df_dalpha == 0, alpha);
alpha_optimal = double(alpha_optimal);

disp('Optimal alpha:');
disp(alpha_optimal);

% Calculate next iteration point using optimal alpha
x_next = double(iter_vector) - alpha_optimal * double(grad_norm);
disp('Next iteration point:');
disp(x_next);

% Calculate function value at next point
f_next = 25*x_next(1)^2 + x_next(2)^2;
disp('Function value at next point:');
disp(f_next);