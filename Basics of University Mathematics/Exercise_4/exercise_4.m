% Exercise 4
% Task 1
clearvars
close all
clc

% a
f = @(x) exp(x);
f(-1)

% b
% syms f(x);
% f(x) = exp(x);
% f(-1)

% c
x = 3;
a = 4;

f(x)
f(a)

% d
h = @(x) x.^2;
x = 1:10;
h(x)

syms h(z);
h(z) = z^2;
h(x)

% e
f = @(x) x.^2;
g = @(x) (x+1).^2-1;

x = linspace(-2, 2, 1000);

fplot(f, [-2, 2]);
figure;
plot(x, g(x));

%% Exercise 4
% Task 2
clearvars
close all
clc

x = 1:1000;
f = @(x) x.^2 + 3*sin(4.*x) - exp(-x);

syms g(z)
g(z) = z^2 + 3*sin(4*z) - exp(-z);

tic
f(x);
numerical_time = toc;

tic
g(x);
symbolic_time = toc;

sym_to_num_ratio = symbolic_time/numerical_time

% b
f = @(x) cos(x) - 1/sqrt(2);

solution1 = fzero(f, pi/4) 
solution2 = fzero(f, 5*pi/4)

%% Exercise 4
% Task 4
clearvars
close all
clc

% a
xi = [2017, 2012, 2008, 2004, 2000];
yi = [7697, 7191, 8952, 8752, 7549];
n = 5;

sum_xiyi = sum(xi.*yi);
sum_xi_sq = sum(xi.^2);

SSxy = sum_xiyi - (sum(xi)*sum(yi)/n)
SSxx = sum_xi_sq - (sum(xi))^2/n

a = SSxy/SSxx
b = mean(yi) - a*mean(xi)

% b
f = @(x) a.*x + b;

plot(xi, yi, '*');
hold on;
plot(xi, f(xi));
hold off;

axis([1995, 2020, 0, 12000]);

% c
title("Number of votes for SDP in the municipal elections in Lappeenranta");
xlabel('$x$-axis (Years)', 'Interpreter', 'latex');
ylabel('$x$-axis (Number of votes)', 'Interpreter', 'latex');
legend({'Scatter plot', 'Fitted line'}, 'Interpreter', 'latex');
grid on;

%% Exercise 4
% Task 5
clearvars
close all
clc

% a
A = [1, 2, 3, 4, 5, 6, 7, 8; 9, 11, 13, 15, 17, 19, 21, 23]

% BONUS
A = [1:8; 9:2:23]

% b
u = A'

% d
A = [3, 2; 3, 5];
B = [4, 3; 1, 2];
x = [4; 1];

A + B'
A*B
B*A
A*x
% x*A 

% The last one is problematic because the columns of the first matrix
% doesn't match the number of rows of the second matrix.

% e
A = [1, 0, 1; 2, -2, 1; 1, 2, 1]
b = [1; 0; 2]
augmented_matrix = [A, b]
rref(augmented_matrix)

x = 0
y = 0.5
z = 1

%% Exercise 4
% Task 6
clearvars
close all
clc

% a i
% By typing "rand" in the command window several times, it seems like the 
% command is returning a random proper fraction value.

% a ii
% rand(m, n) returns an m by n matrix of random proper fraction values.

% a iii A
% Since the default range of values is [0, 1], multiplying it with another
% multiplier 'a' simply changes the range to [0, a]

% a iii B
% Adding a value 'b' with the random number changes the range of values
% from [0, 1] to [b, b+1]

% a iii C
% By the combining both methods of multiplying and adding some constants
% with the rand() function, the range of the random number changes from
% [0,1] to [b, a+b]

% a iv
r1 = rand(1, 10);
r2 = 5 * rand(1, 7);
r3 = 6 + (8 - 6) * rand(1, 15);
r4 = -5 + (5 - (-5)) * rand(15, 10);
r5 = -20 + (-15 - (-20)) * rand(10, 10);
r6 = -0.7 + (0.1 - (-0.7)) * rand(50, 50);

% b
ri = randi([-5, 3], 4, 6);

%% Exercise 4
% Task 7
clearvars
close all
clc

% a
A1 = ones(3, 4); 
A2 = ones(2, 5);
% matrix of any size with all elements set to 1

B1 = zeros(3, 4); 
B2 = zeros(2, 5);
% matrix of any size with all elements set to 0

C1 = eye(3);
C2 = eye(3, 4);
% generates an identity matrix of any size

% if passed only one constant (n) parameter in the function, it generates a 
% square matrix of n by n size.

% b i
V1 = round(10 * rand(1, 5));
V2 = round(-5 + 10 * rand(1, 7));
V3 = round(50 * rand(1, 4));
V4 = round(20 * rand(1, 10));

% b ii
M = [-3.7, 2.5, -1.1, 4.9; 0.3, -5.8, 2.1, -6.5];

R_round = round(M);   % rounds to nearest integer
R_ceil = ceil(M);     % rounds towards positive infinity
R_floor = floor(M);   % rounds towards negative infinity
R_fix = fix(M);       % rounds towards zero

% BONUS
x = 3.14159265359;
rounded_to_3_decimal = round(x, 3);
rounded_to_5_decimal = round(x, 5);
