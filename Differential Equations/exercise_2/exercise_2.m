% Exercise 2
% Task 1
clearvars
clc
close all

% a i
syms x
f(x) = x^2 - 3;
val_at_2 = f(2);
% Using symfun allows direct evaluation f(2), whereas sym would require 
% subs(f, x, 2)

% a ii
sol_eqn = solve(f(x) == 4, x);

% a iii
syms y
g(x,y) = 2*x^2 + 3*y^2;
val_2_3 = g(2,3);
val_x_3 = g(x,3);
sol_y = solve(g(1,y) ==14, y);

% Bonus
% For a i
f_sym = x^2 -3;
bonus_val = subs(f_sym, x, 2);

% For a iii
bonus_expr = subs(g(x,y), x, 1);
bonus_sol = solve(bonus_expr ==14, y);

% b i
syms y(t)
k = -0.2;
ode = diff(y,t) == k*y;
sol_gen = dsolve(ode)

% b ii
cond = y(0) == 2;
sol_ivp = dsolve(ode, cond)

figure;
fplot(sol_ivp, [0, 10]);
hold on;
t_point = 2;
y_point = subs(sol_ivp, t, t_point);    % value at t = 2
plot(t_point, y_point, 'ro');
title('Solution of y'' = -0.2y, y(0)=2');
xlabel('t'); ylabel('y(t)');
legend('Solution', 'Point at t=2');
grid on;
hold off;


%% Exercise 2
% Task 2
clearvars
clc
close all

% a
syms y(t) k C1
ode = diff(y,t) == k*y;
sol_gen = dsolve(ode)

% b
syms y(t) k
cond = y(0) == 2;
ysol(t, k) = dsolve(ode, cond)

% c
syms k
eqn = ysol(30, k) == 1/2;  
k_sol = solve(eqn, k);
disp('2(c) Value of k:');
disp(k_sol);

% d
k_numeric = double(k_sol);
disp('2(d) Numeric k:');
disp(k_numeric);

% e
syms t
eqn_time = ysol(t, k_sol) == 0.02;  % 2*exp(k*t) = 0.02
t_sol = solve(eqn_time, t);
t_numeric = double(t_sol);
disp('2(e) Symbolic time:');
disp(t_sol);
disp('Numeric time (years):');
disp(t_numeric);

%% Exercise 2
% Task 3
clearvars
clc
close all

anim_delay = 0;

% a i
syms x
f = x^2;

% a ii
g = @(x) x.^2;

% a iii
h = matlabFunction(f);

disp(['Type of f: ', class(f)]);  % Symbolic expression
disp(['Type of g: ', class(g)]);  % Function handle
disp(['Type of h: ', class(h)]);  % Function handle

% a iv
x_val = pi;
f_value = subs(f, x, x_val);
h_value = h(x_val);
disp(['Value of f(pi): ', num2str(double(f_value))]);
disp(['Value of h(pi): ', num2str(h_value)]);

% a v
x_vals = 1:10000;

tic;
f_values = subs(f, x, x_vals); 
toc;

tic;
h_values = h(x_vals); 
toc;

% b
tspan = 0:0.01:10;
y0_vals = [-5:1:-1, 1:5]; 

% b i

figure;
hold on;
tic;
for y0 = y0_vals
    [t, y] = ode45(@(t, y) -y, tspan, y0);
    plot(t, y);
    pause(anim_delay);
end
toc;
title('ODE Solution using ode45');
xlabel('t'); ylabel('y(t)');
grid on;
hold off;

% b ii
syms y(t)
y_general = dsolve(diff(y) == -y);
figure;
hold on;
tic;
for y0 = y0_vals
    y_sol = dsolve(diff(y) == -y, y(0) == y0);
    fplot(y_sol, [0, 10]);
    pause(anim_delay);
end
toc;
title('ODE Solution using dsolve');
xlabel('t'); ylabel('y(t)');
grid on;
hold off;

% b iii
syms y(t) C1 t
y_general = dsolve(diff(y) == -y, 'IgnoreAnalyticConstraints', true);

hfun = matlabFunction(y_general, 'Vars', [t, C1]); 

figure;
hold on;
tic;
for y0 = y0_vals
    C1_value = y0; 
    fplot(@(t) hfun(t, C1_value), [0, 10]); 
    pause(anim_delay);
end
toc;
title('ODE Solution using matlabFunction');
xlabel('t'); ylabel('y(t)');
grid on;
hold off;

% b iv
figure;
hold on;
tic;
for y0 = y0_vals
    t_vals = linspace(0, 10, 100);
    y_vals = hfun(t_vals, y0);
    plot(t_vals, y_vals);
    pause(anim_delay);
end
toc;
title('ODE Solution using plot');
xlabel('t'); ylabel('y(t)');
grid on;
hold off;
