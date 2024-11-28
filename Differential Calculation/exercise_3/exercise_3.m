% Exercise 3
% Task 1

% a
x = linspace(-5, 5, 50);
y = linspace(-5, 5, 50);

[X, Y] = meshgrid(x, y);
Z = 2*X + 3*Y + 1;

% Plot using surf
figure;
surf(X, Y, Z);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Surface plot of f(x, y) = 2x + 3y + 1');
colorbar;

% Plot using mesh
figure;
mesh(X, Y, Z);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Mesh plot of f(x, y) = 2x + 3y + 1');
colorbar;

%% b
x = linspace(-3, 3, 100);
y = linspace(-2, 2, 100);

[X, Y] = meshgrid(x, y);

Z = X .* Y .* exp(-X.^2 - Y.^2);

% Plot using surf
figure;
surf(X, Y, Z);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Surface plot of f(x, y) = xye^{-x^2 - y^2}');
colorbar;

% Plot using mesh
figure;
mesh(X, Y, Z);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Mesh plot of f(x, y) = xye^{-x^2 - y^2}');
colorbar;

% c
% surf: creates a surface of the function. The surface is filled.
% mesh: creates a grid of the surface of the function. The grid isn't
%       filled

%% d i
[X, Y] = meshgrid(1:5, 1:3);
% X
% Y
plot(X, Y, 'ro');
hold on;
plot(X(2, 4), Y(2, 4), 'kp');
% X is a matrix that contains the x-coordinates repeated across rows.
% Y is a matrix that contains the y-coordinates repeated across columns.
% First plot command: plots all grid points as red circles.
% Second plot command: Highlights the specific grid point (2, 4) as a black
%                      pentagram.

%% d ii
[X, Y] = meshgrid(1:5, 1:3);
Z = X - Y + 1;
% Z
plot3(X, Y, Z, 'ro');
hold on;
plot3(X(2, 4), Y(2, 4), Z(2, 4), 'kp');
% X and Y are the same as d(i)
% Z computes the z-coordinate for each grid point.
% The first plot3: plots all grid point as red circles in 3D space.
% The second plot3: Hightlights the specific point X(2, 4), Y(2, 4), 
%                   Z(2, 4) as a black pentagram in 3D space.

% d iii
% mesh-command connects all the points in the grid, creating a wireframe
% in the 3D space. 
%% Exercise 3
% Task 2
clearvars
close all
clc

% c
a = input('Enter a number greater than 1: ');
if a <= 1
    error('The number must be greater than 1.');
end

while a <= 1e6
    a = a^2; 
    % fprintf('Current value of a: %.2f\n', a); 
end

fprintf('The final value of a that exceeds one million is %.2f\n', a);

% if a = 1, then the program will fall into an infinite loop.
%% Exercise 3
% Task 3
clearvars
close all
clc

% a
q = input('Enter the common ratio q: ');
n = input('Enter the number of terms n: ');

sum = 0;
i = 0;

while i <= n
    sum = sum + q^i; 
    i = i + 1; 
end

fprintf('The geometric sum is: %.2f\n', sum);

%% b
clearvars;

q = input('Enter the common ratio q: ');
n = input('Enter the number of terms n: ');

terms = q .^ (0:n); 
sum_a = sum(terms);

fprintf('The geometric sum is: %.2f\n', sum_a);

%% Exercise 3
% Task 4
clearvars
close all
clc

% a
data = readmatrix("Differential Calculation\exercise_3\MES3E4.xlsx");
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
sum_of_errors = sum(abs(error))

% iv
% The positive or negative error values might affect the sum. Since we want
% the sum of all the errors, taking the absolutes of each of the error 
% values won't affect the sum in the end.

% e
a = -0.60;
b = -1.36;
x = floor(X(1,1)):ceil(X(end, 1));
y = a*x + b;
plot(x, y, 'g', 'DisplayName', 'Improved Estimated Line');

y_estimated = a*X + b;
error = y_estimated - Y;
sum_of_errors_improved = sum(abs(error))
hold off;
legend('Location', 'best'); 
% Visually, a better fit has been achieved.

%% Exercise 3
% Task 5
clearvars
close all
clc

% c
a = input('Enter a number greater than 1: ');
if a <= 1
    error('The number must be greater than 1.');
end

while true
    a = a^2; 
    if a >= 1e6
        break
    end
end

fprintf('The final value of a that exceeds one million is %.2f\n', a);