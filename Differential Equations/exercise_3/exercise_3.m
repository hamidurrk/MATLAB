%% Exercise 3
% Task 1
clearvars
clc
close all


syms y(t)
Dy = diff(y, t);
D2y = diff(Dy, t);

% b i
eq1 = D2y - y == 0;
conds1 = [y(0) == 1, Dy(0) == 1];
sol1 = dsolve(eq1, conds1);

% b ii
eq2 = D2y + 2*Dy + y == 0;
conds2 = [y(0) == 1, Dy(0) == 1];
sol2 = dsolve(eq2, conds2);

% b iii
eq3 = D2y + 9*y == 0;
conds3 = [y(0) == 1, Dy(0) == 1];
sol3 = dsolve(eq3, conds3);

disp('Solution for (i):'), disp(sol1);
disp('Solution for (ii):'), disp(sol2);
disp('Solution for (iii):'), disp(sol3);

fplot(sol1, [0, 5]); hold on;
fplot(sol2, [0, 5]);
fplot(sol3, [0, 5]);
legend('y'''' - y = 0', 'y'''' + 2y'' + y = 0', 'y'''' + 9y = 0');
title('Solutions to Differential Equations');
xlabel('t');
ylabel('y(t)');
grid on;
hold off;


%% c
clearvars
clc
close all

syms y(t)
Dy = diff(y, t);
D2y = diff(Dy, t);

% iv
eq4 = Dy == sin(y) + sin(t);
sol4 = dsolve(eq4);
disp('Solution for (iv):');
disp(sol4);

% v
eq5 = D2y + sin(t)*Dy + y == 0;
sol5 = dsolve(eq5);
disp('Solution for (v):');
disp(sol5);

% d
% Eqn i:    Linear, second-order, homogenous
% Eqn ii:   Linear, second-order, homogenous
% Eqn iii:  Linear, second-order, homogenous
% Eqn iv:   Non-linear, first-order, non-homogenous
% Eqn v:    Linear, second-order, homogenous

%% Exercise 3
% Task 2
clearvars
clc
close all

% a
syms y(t)
Dy = diff(y, t);
D2y = diff(Dy, t);

eqn = D2y + y == 0;

conds1 = [y(1) == 2, Dy(1) == 0];
conds2 = [y(1) == 0, Dy(1) == 0];
conds3 = [y(1) == -2, Dy(1) == 0];

sol1 = dsolve(eqn, conds1);
sol2 = dsolve(eqn, conds2);
sol3 = dsolve(eqn, conds3);

fplot(sol1, [0, 10], 'r', 'LineWidth', 2); hold on;
fplot(sol2, [0, 10], 'b', 'LineWidth', 2);
fplot(sol3, [0, 10], 'g', 'LineWidth', 2);
legend('y(1) = 2, y''(1) = 0', 'y(1) = 0, y''(1) = 0', 'y(1) = -2, y''(1) = 0');
title('Solutions to y'''' + y = 0 with Different Initial Conditions');
xlabel('t'); ylabel('y(t)');
grid on; hold off;


%% b
clearvars; 
close all;

syms y(t)
Dy = diff(y, t);
D2y = diff(Dy, t);

eqn = D2y + y == 0;

conds1 = [y(1) == 1, Dy(1) == -2];
conds2 = [y(1) == 1, Dy(1) == 0];
conds3 = [y(1) == 1, Dy(1) == 2];

sol1 = dsolve(eqn, conds1);
sol2 = dsolve(eqn, conds2);
sol3 = dsolve(eqn, conds3);

figure;
fplot(sol1, [0, 10], 'r', 'LineWidth', 2); hold on;
fplot(sol2, [0, 10], 'b', 'LineWidth', 2);
fplot(sol3, [0, 10], 'g', 'LineWidth', 2);
legend('y(1) = 1, y''(1) = -2', 'y(1) = 1, y''(1) = 0', 'y(1) = 1, y''(1) = 2');
title('Effect of y''(a) on Solution');
xlabel('t'); ylabel('y(t)');
grid on; hold off;


%% c
clearvars; 
close all;


syms y(t)
Dy = diff(y, t);
D2y = diff(Dy, t);

eqn = D2y + y == 0;

syms C1 C2
y_sol = C1*cos(t) + C2*sin(t);

eqs1 = subs(y_sol, t, 1) == 1;
eqs2 = subs(y_sol, t, 10) == -2;
sols = solve([eqs1, eqs2], [C1, C2]);
sol1 = subs(y_sol, [C1, C2], [sols.C1, sols.C2]);

eqs3 = subs(y_sol, t, 1) == -3;
eqs4 = subs(y_sol, t, 10) == 2;
sols2 = solve([eqs3, eqs4], [C1, C2]);
sol2 = subs(y_sol, [C1, C2], [sols2.C1, sols2.C2]);

figure;
fplot(sol1, [0, 10], 'r', 'LineWidth', 2); hold on;
fplot(sol2, [0, 10], 'b', 'LineWidth', 2);
legend('y(1) = 1, y(10) = -2', 'y(1) = -3, y(10) = 2');
title('Solutions to y'''' + y = 0 with Boundary Conditions');
xlabel('t'); ylabel('y(t)');
grid on; hold off;


%% Exercise 3
% Task 3
clearvars
clc
close all


% a
x = linspace(0, 2*pi, 100); 
y = sin(x);

figure;
hold on;
plot(x, y, 'b', 'LineWidth', 1.5); % Plot the sine function
h = plot(x(1), y(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Marker

for i = 1:length(x)
    set(h, 'XData', x(i), 'YData', y(i));
    pause(0.05); % Pause for animation effect
end

hold off;
%% b
clc; clear; close all;
syms s(t) k m g

m_val = 1;       % Mass (kg)
k_val = 0.5;     % Air resistance coefficient
g_val = 9.81;    % Gravitational acceleration (m/s²)
s0 = 100;        % Initial position (m)
v0 = 0;          % Initial velocity (m/s)
v_up = 20;       % Initial velocity (m/s) (upward throw)
v_down = -20;    % Initial velocity (m/s) (downward throw)

eq1 = m * diff(s, t, 2) + m * g == 0;  % Without air resistance
eq2 = m * diff(s, t, 2) + k * diff(s, t) + m * g == 0; % With air resistance

Ds = diff(s, t);

sol1 = dsolve(eq1, s(0) == s0, Ds(0) == v0);  
sol2 = dsolve(eq2, s(0) == s0, Ds(0) == v0);  
sol3 = dsolve(eq2, s(0) == s0, Ds(0) == v_up);
sol4 = dsolve(eq2, s(0) == s0, Ds(0) == v_down);

sol1_func = matlabFunction(subs(sol1, [m, g], [m_val, g_val]), 'Vars', t);
sol2_func = matlabFunction(subs(sol2, [k, m, g], [k_val, m_val, g_val]), 'Vars', t);
sol3_func = matlabFunction(subs(sol3, [k, m, g], [k_val, m_val, g_val]), 'Vars', t);
sol4_func = matlabFunction(subs(sol4, [k, m, g], [k_val, m_val, g_val]), 'Vars', t);

t_vals = linspace(0, 5, 100);
s1_vals = sol1_func(t_vals);
s2_vals = sol2_func(t_vals);
s3_vals = sol3_func(t_vals);
s4_vals = sol4_func(t_vals);

figure;
hold on;
plot(t_vals, s1_vals, 'r', 'LineWidth', 1.5); % No Air Resistance
plot(t_vals, s2_vals, 'b', 'LineWidth', 1.5); % With Air Resistance
plot(t_vals, s3_vals, 'g', 'LineWidth', 1.5); % Thrown Upward
plot(t_vals, s4_vals, 'm', 'LineWidth', 1.5); % Thrown Downward

h1 = plot(t_vals(1), s1_vals(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
h2 = plot(t_vals(1), s2_vals(1), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
h3 = plot(t_vals(1), s3_vals(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
h4 = plot(t_vals(1), s4_vals(1), 'mo', 'MarkerSize', 10, 'MarkerFaceColor', 'm');

xlabel('Time (s)');
ylabel('Position (m)');
legend('No Air Resistance', 'With Air Resistance', 'Thrown Upward', 'Thrown Downward');
grid on;
hold off;

% Animation Loop
for i = 1:length(t_vals)
    set(h1, 'XData', t_vals(i), 'YData', s1_vals(i));
    set(h2, 'XData', t_vals(i), 'YData', s2_vals(i));
    set(h3, 'XData', t_vals(i), 'YData', s3_vals(i));
    set(h4, 'XData', t_vals(i), 'YData', s4_vals(i));
    
    pause(0.05);
end

