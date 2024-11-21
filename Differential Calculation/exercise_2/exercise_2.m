% Exercise 2
% Task 1

% b
a = 0;
b = 1;
if a == b
    c = 1;
elseif a > b + 2
    c = 2;
elseif a > b
    c = 3;
else
    c = 4;
end
% for 1, 4, 2, and 0 values of the variable 'a', 
% the value of 'c' would be 1, 2, 3, 4 respectively.
c

%% Exercise 2
% Task 2
clearvars
close all
clc

% a
a = input("Please enter a positive integer: ");
a_str = num2str(a);

if a < 0
    disp(['The given number ', a_str, ' is not positive.']);
else
    if isprime(a)
        disp(['The given number ', a_str, ' is a prime.']);
    else
        disp(['The given number ', a_str, ' is not a prime.'])
    end
end

% b
a = input("Please enter a positive integer: ");
a_str = num2str(a);

if a < 0
    disp(['The given number ', a_str, ' is not positive.']);
elseif a > 0 && isprime(a)
    disp(['The given number ', a_str, ' is a prime.']);
else
    disp(['The given number ', a_str, ' is not a prime.'])
end

%% Exercise 2
% Task 3
clearvars
close all
clc

weight = input('Enter your weight in kilograms (kg): ');
height = input('Enter your height in meters (m): ');

if weight <= 0 || height <= 0
    disp('Error: Please check the entered values. Weight and height must be positive numbers.');
else
    bmi = weight / (height^2);
    
    fprintf('Your BMI is: %.2f\n', bmi);
    
    if bmi < 18.5
        disp('Category: Underweight');
    elseif bmi >= 18.5 && bmi < 24.9
        disp('Category: Normal weight');
    elseif bmi >= 25 && bmi < 29.9
        disp('Category: Overweight');
    else
        disp('Category: Obese');
    end
end

%% Exercise 2
% Task 5
clearvars
close all
clc

% a (numeric)
n = @(x) cos(x);
values_numeric = n([0, pi/4, pi])

% a (symbolic)
syms x
s = cos(x);
values_symbolic = subs(s, x, [0, pi/4, pi])

% b
x_numeric = linspace(-6, 6, 1000);
y_numeric = n(x_numeric);

figure;
plot(x_numeric, y_numeric, 'b--', 'LineWidth', 2);
hold on;

fplot(n, [-6, 6], 'r--');
hold off;

xlabel('x');
ylabel('cos(x)');
legend('Numeric (plot)', 'Symbolic (fplot)');
title('Comparison of Numeric and Symbolic Plots');
grid on;

% c
syms x h

limit1 = limit((cos(x) - 1)/x, x, 0);                       % i
limit2 = limit((log(x + h) - log(x))/h, h, 0);              % ii
limit3 = limit((2*x - 6*x^4)/(-7*x^3 + x^2), x, inf);       % iii
limit4 = limit((x^2 - 2*x - 15)/(x^2 - 10*x + 25), x, 5);   % iv
% the result of part (ii) was 1/x

% BONUS
bonus_limit1 = limit(abs(16 - x)/(4 - sqrt(x)), x, 16, 'right');
bonus_limit2 = limit(abs(16 - x)/(4 - sqrt(x)), x, 16, 'left');

%% Exercise 2
% Task 6
clearvars
close all
clc

% a
syms y(x)

dy_dx = diff(y(x));
% Result: diff(y(x), x)

d_y_squared_dx = diff(y(x)^2);
% Result: 2*y(x)*diff(y(x), x)

% In standard mathematical notation:
% dy/dx
% d(y^2)/dx = 2*y*dy/dx

% b
syms y(x)
eq = x^2 + y^2 == 4;

d_eq = diff(eq, x);

dy_dx = solve(d_eq, diff(y, x));

% dy/dx = -x/y

% c
syms x(t) y(t)
eq = x(t)^2 + y(t)^2 == 4;

d_eq = diff(eq, t); 

dy_dt = solve(d_eq, diff(y(t), t));

% dy/dt = -(x*dx/dt) / y

%% Exercise 2
% Task 7
clearvars
close all
clc

% a
t = linspace(0, 40, 1000); 

x = (t / 10) .* cos(t);
y = (t / 10) .* sin(t);

figure;
plot(x, y, 'b-', 'LineWidth', 1.5);
xlabel('x(t)');
ylabel('y(t)');
title('Parametric Curve: x(t) = (t/10)cos(t), y(t) = (t/10)sin(t)');
grid on;
axis equal;

% b
t = linspace(0, 2*pi, 1000); 

x = cos(t);
y = sin(t);

figure;
plot(x, y, 'r-', 'LineWidth', 1.5);
xlabel('x(t)');
ylabel('y(t)');
title('Parametric Curve: f(t) = cos(t)i + sin(t)j');
grid on;
axis equal;

% c
t = linspace(0, 40, 1000);

x = (t / 10) .* cos(t);
y = (t / 10) .* sin(t);
z = t / 10;

figure;
plot3(x, y, z, 'g-', 'LineWidth', 1.5);
xlabel('x(t)');
ylabel('y(t)');
zlabel('z(t)');
title('3D Parametric Curve: f(t) = (t/10)cos(t)i + (t/10)sin(t)j + (t/10)k');
grid on;
axis equal;