% Exercise 3
% Task 1
clearvars
close all
clc

% b
v = [1 2, 3];

% c
x1 = 1:4;

% d
v = [1; 2; 3; 4];
v = transpose(x1);
v = x1';

% e
v1 = 1:10;
v2 = [1:10];
v3 = 1:2:9;
v4 = 1:2:10;
v5 = -25:3:0;
v6 = 1000:-10:950;
v7 = [1:5, 5:-1:0];

vi      = 5:5:125;
vii     = 125:-5:5;
viii    = -100:1:-51;
viv     = 3.1:0.2:43.1;
vv      = [0:100, 100:-1:0, 1:49];
vvi     = [0:2:62, 60:-2:-4];

% f
f1 = [sum(vi), sum(vii), sum(viii), sum(viv), sum(vv), sum(vvi)];
f2 = mean(f1)

% g
g1 = linspace(0, 100, 51);
% There are 51 elements in the vector. The interval [0, 100] is divided 
% into 50 subintervals, because with 51 points, there are 50 spaces 
% between them.

g2 = linspace(0, 100);
% With 2 parameters, the linspace command divides the interval into 100
% points by default.

g3 = linspace(3, 48, 118);
% Since we want 117 (n) equal parts, we would need 118 (n+1) equally
% distributed points to be created. And these 118 points are the endpoints
% of the subintervals that we want.

%% Exercise 3
% Task 2
clearvars
close all
clc

% b
x1 = linspace(0, 2*pi, 100);
sin(x1);

% c
sin(pi:pi:5*pi);
% sin pi to 5pi with an increment of pi should be zero but since the pi in
% MATLAB in just an approxiamation of the actual pi, the output value tends
% to zero but not zero. All of the values are multiplied with 1x10^-15,
% that's why.

% d
y1 = x1.*x1;
% In MATLAB, the * operator is for matrix multiplication, which means it 
% expects the number of columns in the first vector to match the number 
% of rows in the second vector. Trying to multiply x1*x1 directly is an 
% not valid because the dimensions don’t align for matrix multiplication.

% However, the .* operator is for element-wise multiplication. This 
% operator multiplies each element of x1 with the corresponding element 
% of x1. So, x1.*x1 works as it performs element-wise multiplication.

% e
x2 = [0, 1, 2, 3];
x3 = [1, 2, 3, 4];

x2 + 1;
% x2*x3;
x2.*x3;
% x2^2;
x2.^2;
% 2^x2;
2.^x2;
12./x3;
x2./x3;
x2 - x3;
% x2 .+ x3;
2*x2;
2.*x2;

% f
x = 1:10;
y = 1./x;

%% Exercise 3
% Task 3
clearvars
close all
clc

x = [-pi/3, -pi/4, 0, pi/4, pi/3];

% i
lhs1 = sin(3*x);
rhs1 = 3*sin(x) - 4*sin(x).^3;
diff1 = lhs1 - rhs1

% ii
lhs2 = abs(sin(x/2));
rhs2 = sqrt((1 - cos(x))/2);
diff2 = lhs2 - rhs2

% iii
lhs3 = tan(3*x);
rhs3 = (3*tan(x) - tan(x).^3) ./ (1 - 3*tan(x).^2);
diff3 = lhs3 - rhs3

% (BONUS)
a = [-3, -1, 0, 1, 3];
b = [-2, -1, 0, -1, 2];
lhs_bonus = abs(a + b);
rhs_bonus = abs(a) + abs(b);
diff_bonus = lhs_bonus - rhs_bonus
% Here, the equation is not always true. The correct relation would be 
% that: abs(a + b) <= abs(a) + abs(b);

%% Exercise 3
% Task 3
clearvars
close all
clc

% b
x1 = [-3,-1, 0, 1, 3];
y1 = [7, 3, 7, 3, 7];

% plot(x1, y1, '*');
% figure;
% plot(x1, y1, 'o')
% figure
% plot(x1, y1)
% figure
% plot(x1, y1, 'r')
% figure
% plot(x1, y1, 'b')
% figure
% plot(x1, y1, '--')
% figure
% plot(x1, y1, 'ro')
% x2 =-5:0.5:5;
% y2 = sin(x2);
% figure
% plot(x2, y2, 'o')
% figure
% plot(x2, y2)

% b(i)
x1 = [0 1 2];
y1 = [5 4 3];

figure;
plot(x1, y1, 'ro');
figure;
plot(x1, y1, 'r-o');  

% b(ii)
x2 = [-2 5];
y2 = [3 6];

figure;
plot(x2, y2, 'b-');

% b(iii)
x3 = -5:0.01:5;
y3 = sin(x3);

figure;
plot(x3, y3);

%% Exercise 3
% Task 4
clearvars
close all
clc

% a
x1 = -3.2:0.8:3.2;
x2 = -3.2:0.1:3.2;

f1 = x1.^2;
f2 = x2.^2;

figure;
plot(x1, f1, 'ro-');
hold on;
plot(x2, f2, 'b-');
hold off;

title('$f(x) = x^2$ with Different Spacing', 'Interpreter', 'latex');
xlabel('$x$-axis', 'Interpreter', 'latex');
ylabel('$f(x) = x^2$', 'Interpreter', 'latex');
legend({'Points spaced 0.8 units', 'Points spaced 0.1 units'}, 'Interpreter', 'latex');
grid on;

% b
x = -2*pi:0.01:2*pi;

f = sin(x);
g = sin(2*x);
h = sin(3*x);

figure;
plot(x, f, 'r-');
hold on;
plot(x, g, 'g--');
plot(x, h, 'b-.');
hold off;

title('Sine Functions: $f(x) = \sin(x)$, $g(x) = \sin(2x)$, $h(x) = \sin(3x)$', 'Interpreter', 'latex');
xlabel('X-axis', 'Interpreter', 'latex');
ylabel('Y-axis', 'Interpreter', 'latex');
legend({'$f(x) = \sin(x)$', '$g(x) = \sin(2x)$', '$h(x) = \sin(3x)$'}, 'Interpreter', 'latex');
grid on;

% BONUS - used LaTeX

%% Exercise 3
% Task 5
clearvars
close all
clc

u = [1, 2, 3];
v = [-2, -3, 1];

% a
sum_uv = u + v
diff_uv = u - v
linear_combination = 2*u - 3*v

% b
magnitude_u = norm(u)
magnitude_v = norm(v)

% c
unit_u = u / magnitude_u
unit_v = v / magnitude_v

% d
dot_product = dot(u, v);
angle_rad = acos(dot_product / (magnitude_u * magnitude_v));
angle_deg = rad2deg(angle_rad)

% BONUS
function result = customDot(a, b)
    result = sum(a .* b);
end

dot_product_custom = customDot(u, v)

%% Exercise 3
% Task 6
clearvars
close all
clc

u = [1, 2, 3];
v = [-2, -3, 1];

% a
quiver3(0, 0, 0, 1, 1, 1, 'AutoScale', 'off');

% b
cross_uv = cross(u, v);
cross_vu = cross(v, u);

origin = [0, 0, 0];

figure;
hold on;

quiver3(origin(1), origin(2), origin(3), u(1), u(2), u(3), 'r', 'AutoScale', 'off');
quiver3(origin(1), origin(2), origin(3), v(1), v(2), v(3), 'g', 'AutoScale', 'off');
quiver3(origin(1), origin(2), origin(3), cross_uv(1), cross_uv(2), cross_uv(3), 'b', 'AutoScale', 'off');
quiver3(origin(1), origin(2), origin(3), cross_vu(1), cross_vu(2), cross_vu(3), 'm', 'AutoScale', 'off');

axis([-5 5 -5 5 -5 5]); 

view(3);

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Vector Plot');
grid on; 

% BONUS
plot3([-5, 5], [0, 0], [0, 0], 'k');
plot3([0, 0], [-5, 5], [0, 0], 'k');
plot3([0, 0], [0, 0], [-5, 5], 'k');

hold off;

%% Exercise 3
% Task 7
clearvars
close all
clc

% a
% f(x) = x^2
% fplot(f(x))

% F9 (executes the selected line) and F5 (runs the entire script)
% The function here isn't declared properly. The correct syntax has been
% shown below
% 
f = @(x) x.^2
fplot(f)


% b
% x = [1, 2, 2, 3, 1.5, 0, 1];
% y = [0, 0, 3, 3, 4, 3, 3, 0];
% plot(x, y)

% The error is due to the y vector having an incorrect number of elements 
% compared to the x vector

x = [1, 2, 2, 3, 1.5, 0, 1, 1];
y = [0, 0, 3, 3, 4, 3, 3, 0];
figure
plot(x, y, 'k');


% c
x = -5:0.1:5;
y = x.^2;
plot(y, 'k')

% To plot y with respect to x we will have to write the code accordingly.
% The plot command is missing an argument here.
% 
% x = -5:0.1:5; 
% y = x.^2;
% figure
% plot(x, y, 'k');
