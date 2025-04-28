% Constants and initial values
A = 3; 
B = -2; 
C = 3;

x_0 = 5; 
y_0 = 7;

iter_vector = [x_0; y_0]; 

% Function
syms x y alpha;
f = A*x^2 - B*x*y + C*y^2 + x - y; 
grad_f = gradient(f, [x, y]); 

grad_val = subs(grad_f, [x, y], iter_vector'); 
grad_norm = grad_val / norm(grad_val);

x_next = iter_vector - alpha * grad_norm;

% f(x_next)
f_alpha = subs(f, [x, y], x_next');

% Dynamic alpha calculation
df_dalpha = diff(f_alpha, alpha);
alpha_optimal = solve(df_dalpha == 0, alpha);

% Calculate next iteration point
x_next = iter_vector - alpha_optimal * grad_norm;
f_next = subs(f, [x, y], x_next');

