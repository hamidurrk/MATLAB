%% MATLAB 1
clear; clc; close all;

% Coefficients
a = 1;
b = 3;
c = -16;

% 1. Discriminant
D = b^2 - a*c;
disp(['Discriminant D = ', num2str(D)])

if D > 0
    disp('The PDE is HYPERBOLIC');
end

% 2. Characteristic slopes
p_plus  = (b + sqrt(D)) / a;
p_minus = (b - sqrt(D)) / a;

disp(['p+ = ', num2str(p_plus)])
disp(['p- = ', num2str(p_minus)])

% 3. Plot characteristic line families
x = linspace(-5, 5, 400);

figure; hold on; grid on;
title('MATLAB 1: Characteristic Families')

Cvals = -5:1:5;

for C = Cvals
    y1 = p_plus*x + C;
    y2 = p_minus*x + C;

    plot(x, y1, 'b', 'LineWidth', 1.5)
    plot(x, y2, 'r--', 'LineWidth', 1.5)
end

xlabel('x'); ylabel('y');
legend('y - p_+ x = C_+', 'y - p_- x = C_-');

% 4. Characteristic coordinates and level sets
[X, Y] = meshgrid(linspace(-5,5,200), linspace(-5,5,200));

xi  = Y - p_plus  .* X;
eta = Y - p_minus .* X;

figure; hold on; grid on;
title('Level Sets of Characteristic Coordinates')

contour(X, Y, xi, 15, 'b')
contour(X, Y, eta, 15, 'r--')

xlabel('x'); ylabel('y');
legend('\xi = const', '\eta = const');


%% MATLAB 2 
clear; clc; close all;

% Coefficients
a = 1;
b = -3;
c = 9;

% 1. Discriminant
D = b^2 - a*c;
disp(['Discriminant D = ', num2str(D)])

if D == 0
    disp('The PDE is PARABOLIC');
end

% 2. Repeated characteristic slope
p = b / a;
disp(['Repeated slope p = ', num2str(p)])

x = linspace(-5, 5, 400);

figure; hold on; grid on;
title('MATLAB 2: Single Characteristic Family')

Cvals = -5:1:5;

for C = Cvals
    y = p*x + C;
    plot(x, y, 'b', 'LineWidth', 1.5)
end

xlabel('x'); ylabel('y');
legend('y - px = C');

% 3. Characteristic coordinates and Jacobian
[X, Y] = meshgrid(linspace(-5,5,200), linspace(-5,5,200));

eta = Y - p .* X;   % characteristic coordinate
xi  = X;            % second coordinate

% Jacobian determinant
% xi_x = 1, xi_y = 0
% eta_x = -p, eta_y = 1
J = 1*1 - 0*(-p);
disp(['Jacobian J = ', num2str(J)])

% Plot level sets
figure; hold on; grid on;
title('Level Sets in Parabolic Case')

contour(X, Y, eta, 15, 'b')  % characteristic family
contour(X, Y, xi, 15, 'r--') % vertical lines

xlabel('x'); ylabel('y');
legend('\eta = const (characteristic)', '\xi = const');
