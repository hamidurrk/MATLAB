% Exercise 1
% Task 1
clearvars
clc
close all

% a
odefun = @(t, y) - 0.2*y;
tspan = [0, 5];
y0 = 1.5;

[T, Y] = ode45(odefun, tspan, y0); 

% odefun: Defines the ODE y′ = −0.2y using an anonymous function.
% tspan: Specifies the interval [0,5] over which the solution is computed.
% y0: Sets the initial condition y(0) = 1.5
% [T, Y]: Stores the output, where T contains the time points and Y the 
%           corresponding numerical solution for y(t).

% b
figure;
plot(T, Y, 'b-', 'LineWidth', 1.5);
title('Numerical Solution of y'' = -0.2y');
grid on;

% c
y_analytical = @(t) 1.5 * exp(-0.2 * t);
Y_analytical = y_analytical(T);

hold on;
plot(T, Y_analytical, 'r.', 'LineWidth', 1.5); 
legend('Numerical Solution (ode45)', 'Analytical Solution', 'Location', 'northeast');
hold off;

%% d
clearvars
clc
close all

odefun = @(t, y) 0.4 * y; % rewritten as y' = 0.4y
tspan = [0, 3];
y0 = 0.6;

[T, Y] = ode45(odefun, tspan, y0);

% analytical solution
y_analytical = @(t) 0.6 * exp(0.4 * t);
Y_analytical = y_analytical(T);

figure;
plot(T, Y, 'b-', 'LineWidth', 1.5); 
hold on;
plot(T, Y_analytical, 'r.', 'LineWidth', 1.5); 
xlabel('Time t');
ylabel('Solution y(t)');
title('Solution of y'' - 0.4y = 0');
legend('Numerical Solution (ode45)', 'Analytical Solution', 'Location', 'northwest');
grid on;
hold off;

%% e
clearvars
clc
close all

odefun = @(t, y) 0.4 * y;
tspan = [-3, 3];
y0 = 0.6;

[T, Y] = ode45(odefun, tspan, y0);

% analytical solution
y_analytical = @(t) 0.6 * exp(0.4 * t + 1.2);
Y_analytical = y_analytical(T);

figure;
plot(T, Y, 'b-', 'LineWidth', 1.5); 
hold on;
plot(T, Y_analytical, 'r.', 'LineWidth', 1.5); 
xlabel('Time t');
ylabel('Solution y(t)');
title('Solution of y'' - 0.4y = 0 over [-3, 3]');
legend('Numerical Solution (ode45)', 'Analytical Solution', 'Location', 'northwest');
grid on;
hold off;

%% Exercise 1
% Task 2
clearvars
close all
clc

% a
odefun_a = @(x, y) -x * y;
xspan_a = [0, 5];
y0_a = 1; 
[X_a, Y_a] = ode45(odefun_a, xspan_a, y0_a);

y_analytical_a = @(x) y0_a * exp(-x.^2 / 2);
Y_analytical_a = y_analytical_a(X_a);

figure;
plot(X_a, Y_a, 'b-', 'LineWidth', 1.5); 
hold on;
plot(X_a, Y_analytical_a, 'r.', 'LineWidth', 1.5); 
xlabel('x');
ylabel('y(x)');
title('Solution for y'' + xy = 0');
legend('Numerical Solution (ode45)', 'Analytical Solution', 'Location', 'northeast');
grid on;
hold off;

% b
odefun_b = @(x, y) 1 - y;
xspan_b = [0, 5];
y0_b = 2; 
[X_b, Y_b] = ode45(odefun_b, xspan_b, y0_b);

C_b = y0_b - 1; 
y_analytical_b = @(x) 1 + C_b * exp(-x);
Y_analytical_b = y_analytical_b(X_b);

figure;
plot(X_b, Y_b, 'b-', 'LineWidth', 1.5); 
hold on;
plot(X_b, Y_analytical_b, 'r.', 'LineWidth', 1.5); 
xlabel('x');
ylabel('y(x)');
title('Solution for y'' + y = 1');
legend('Numerical Solution (ode45)', 'Analytical Solution', 'Location', 'northeast');
grid on;
hold off;

% c
odefun_c = @(x, y) x - y;
xspan_c = [0, 5];
y0_c = 0; 
[X_c, Y_c] = ode45(odefun_c, xspan_c, y0_c);

C_c = (y0_c + 1) * exp(0); 
y_analytical_c = @(x) x - 1 + C_c ./ exp(x);
Y_analytical_c = y_analytical_c(X_c);

figure;
plot(X_c, Y_c, 'b-', 'LineWidth', 1.5); 
hold on;
plot(X_c, Y_analytical_c, 'r.', 'LineWidth', 1.5); 
xlabel('x');
ylabel('y(x)');
title('Solution for y'' + y = x');
legend('Numerical Solution (ode45)', 'Analytical Solution', 'Location', 'northwest');
grid on;
hold off;

%% Exercise 1
% Task 3
clearvars
close all
clc

% a
f = @(t, y) 0.4 * y; 
y0 = 0.6; 
tspan_a = [-2, 1]; 

[t_a, y_a] = ode45(f, tspan_a, y0);

t_exact = linspace(-2, 1, 100);
y_exact = 0.6 * exp(0.4 * (t_exact + 2));

figure;
plot(t_a, y_a, 'b', 'LineWidth', 1.5); hold on;
plot(t_exact, y_exact, 'r.', 'LineWidth', 1.5);
xlabel('t');
ylabel('y(t)');
legend('Numerical Solution', 'Analytical Solution');
title('Solution of y'' = -0.4y');
grid on;

% b
% It is not possible to solve the problem numerically over [-10, -5]
% because the initial condition is given at t = -2, outside this range.

% c
tspan_c = -2:0.1:1; 
[t_c, y_c] = ode45(f, tspan_c, y0);

if ~ismember(0.9, t_c)
    tspan_c = [-2:0.1:0.9, 0.9:0.1:1]; 
    [t_c, y_c] = ode45(f, tspan_c, y0);
end

% d
index_t1 = find(t_c == 1);
y_at_t1 = y_c(index_t1);

figure;
plot(t_c, y_c, 'b', 'LineWidth', 1.5); hold on;
plot(1, y_at_t1, 'r.', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('t');
ylabel('y(t)');
legend('Numerical Solution', 'Point at t = 1');
title('Solution of y'' = -0.4y with t = 1');
grid on;

% e
index_t09 = find(abs(t_c - 0.9) < 1e-10);
y_at_t09 = y_c(index_t09);

fprintf('Value of y at t = 1: %.4f\n', y_at_t1);
fprintf('Value of y at t = 0.9: %.4f\n', y_at_t09);
