%% MATLAB 1

clear; clc;

% Domain
x = linspace(0, 2, 200);
y = linspace(-6, 6, 200);
[X, Y] = meshgrid(x, y);

% (a) 
% Characteristics: y - 3x = constant
% u is constant along characteristics
U = cos(Y - 3*X);

% (b)
figure;
surf(X, Y, U);
shading interp;
colormap turbo;
hold on;

xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Solution of 2u_x + 6u_y = 0');

% (c)
Cvals = [-5, 0, 5];   % constants for characteristic curves

for C = Cvals
    x_line = linspace(0, 2, 200);
    y_line = 3*x_line + C;
    u_line = cos(y_line - 3*x_line);
    plot3(x_line, y_line, u_line, 'k', 'LineWidth', 2);
end

hold off;

%% MATLAB 2
clc; clear; close all;

% (a) Exact solution from characteristics
% u(x,y) = sin((x - y^2/2)/5) + y/3

% (b) 

% Domain
x = linspace(-10,10,200);
y = linspace(0,5,200);
[X,Y] = meshgrid(x,y);

% Solution
U = sin((X - 0.5*Y.^2)/5) + Y/3;

figure
surf(X,Y,U)
shading interp
hold on
xlabel('x')
ylabel('y')
zlabel('u(x,y)')
title('Solution surface and characteristic curve')
colorbar

% (c) 

% u(0,0) = 0
s = linspace(0,5,200);

x_c = 0.5*s.^2;   % x(s)
y_c = s;          % y(s)
u_c = s/3;        % u(s)

plot3(x_c,y_c,u_c,'r','LineWidth',3)

legend('u(x,y)','Characteristic through (0,0,u(0,0))')
hold off
