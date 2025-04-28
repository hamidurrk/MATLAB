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

alpha = 0.5;

iter_result = double(iter_vector - alpha * grad_norm)
