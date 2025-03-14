% Exercise 5
% Task 1

% b
q = input('Enter the common ratio q: ');
n = input('Enter the number of terms n: ');
v = 1:n;

sum = 0;

for i = 0:n
    sum = sum + q^i; 
    v(i+1) = sum;
end

fprintf('The geometric sum of first n+1 terms is: %.2f\n', sum);
disp(v);

%% Exercise 5
% Task 2
clearvars
close all
clc

% a
n = input('Enter the number of terms in Fibonacci sequece: ');
F = [1, 1];
for i = 3:n
    F(i) = F(i-1)+F(i-2);
end

% b
ratio = 1:n-1;
for i = 1:n-1
    ratio(i) = F(i+1)/F(i);
end
disp(ratio);

% c
golden_ratio = (1+sqrt(5))/2;
diff = [abs(ratio-golden_ratio)];

disp("The absolute values of the differences between the ratios from part (b) and the golden ratio")
disp(diff)

%% Exercise 5
% Task 3
clearvars
close all
clc

% b
theta = linspace(0, 2*pi, 100); 
x = cos(theta); 
y = sin(theta); 

figure; 
plot(x, y, 'b');
hold on; 
axis equal; 
title('Unit Circle with Animated Object');
xlabel('x');
ylabel('y');

% c
t0 = 0; 
x0 = cos(t0); 
y0 = sin(t0); 
% plot(x0, y0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
resolution = pi/100;

% d
for t = 0:200 
    cla; 
    plot(x, y, 'b'); 
    xMarker = cos(t*resolution); 
    yMarker = sin(t*resolution);
    plot(xMarker, yMarker, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
    pause(0.05); 
end
hold off; 



%% Exercise 5
% Task 4
clearvars
close all
clc

% a
theta = linspace(0, 2*pi, 100); 
x = cos(theta); 
y = sin(theta); 

figure; 
plot(x, y, 'b');
hold on; 
axis equal; 
title('Unit Circle with Animated Object');
xlabel('x');
ylabel('y');

% b
t0 = 0; 
x0 = cos(t0); 
y0 = sin(t0); 
h = plot(x0, y0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
resolution = pi/100;

% c, d
for t = 0:200 
    xMarker = cos(t*resolution);
    yMarker = sin(t*resolution);
    h.XData = xMarker;
    h.YData = yMarker;
    pause(0.01); 
end

%% Exercise 4
% Task 5
clearvars
close all
clc

% a
[x, y] = meshgrid(-3:0.1:3, -3:0.1:3);
f = x - y;

figure; 
surf(x, y, f); 
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
title('Surface Plot of f(x, y) = x - y');

% b
figure; 
surf(x, y, f, 'FaceAlpha', 0.1); 
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
title('Translucent Surface Plot of f(x, y) = x - y');

% c
level_curves = [0, 1, 2, 3, -1, -2, -3]; 
figure; 
hold on; 

for k = 1:length(level_curves)
    contour(x, y, f, [level_curves(k), level_curves(k)], 'LineWidth', 2); 
    pause(0.5);
end

% d
axis([-3, 3, -3, 3]);

% e
% The level curve corresponding to the zeros of the function is the one where f(x, y) = 0.
% This occurs when x - y = 0, which is the line y = x in the xy-plane.


%% Exercise 4
% Task 6
clearvars
close all
clc

% a
data = readmatrix("C:\Users\hamid\OneDrive\Documents\MATLAB\Differential Calculation\exercise_3\MES3E4.xlsx");
X = data(:, 1);
Y = data(:, 2);

% b
figure;
plot(X, Y, 'ro', 'DisplayName', 'Original Data');
grid on;
hold on;

a = -0.6;
b = -1.3;
x = floor(X(1,1)):ceil(X(end, 1));
y = a*x + b;
plot(x, y, 'b--', 'DisplayName', 'Initial Estimated Line');

xlabel('X-axis Label');
ylabel('Y-axis Label');
title('Linear Fit and Error Analysis');
legend;

% d
% i
y_estimated = a*X + b;
plot(X, y_estimated, 'bo', 'DisplayName', 'Estimated Data Points');

% ii
error = y_estimated - Y;
% These values represent the distances between each point in the plot.

% iii
sum_of_squares = sum(error.^2);

% iv
% The positive or negative error values might affect the sum. Since we want
% the sum of all the errors, taking the absolutes of each of the error 
% values won't affect the sum in the end.

% e
p = polyfit(X, Y, 1);

x = floor(X(1,1)):ceil(X(end, 1));
y = polyval(p, x);

plot(x, y, 'g', 'DisplayName', 'Improved Estimated Line');

y_estimated = a*X + b;
error = y_estimated - Y;
sum_of_errors_improved = sum(abs(error))
hold off;
legend('Location', 'best'); 
% Visually, a better fit has been achieved.

%% Exercise 4
% Task 7
clearvars
close all
clc

% a
data = readmatrix("Differential Calculation\exercise_3\MES3E4.xlsx");
X = data(:, 1);
Y = data(:, 2);

syms a b

% b
y_estimated = a*X + b; 
error = y_estimated - Y; 
S = sum(error.^2); 

% c
grad_S = gradient(S, [a, b]);

% d
solution = solve(grad_S == 0, [a, b]);
a_opt = double(solution.a);
b_opt = double(solution.b);

% e
fprintf('Optimal Slope (a): %.4f\n', a_opt);
fprintf('Optimal Intercept (b): %.4f\n', b_opt);

figure;
plot(X, Y, 'ro', 'DisplayName', 'Original Data');
hold on;
grid on;

x_range = floor(min(X)):ceil(max(X));
y_initial = -0.6 * x_range - 1.3;
plot(x_range, y_initial, 'b--', 'DisplayName', 'Initial Estimated Line');

y_optimized = a_opt * x_range + b_opt;
plot(x_range, y_optimized, 'g', 'DisplayName', 'Optimized Line');

xlabel('X-axis Label');
ylabel('Y-axis Label');
title('Optimized Linear Fit');
legend('Location', 'best');
hold off;
