% Exercise 5
% Task 1
clearvars
close all
clc

% b
rng(1)
r1 = rand(1, 10);
r2 = 5 * rand(1, 7);
r3 = 6 + (8 - 6) * rand(1, 15);
r4 = -5 + (5 - (-5)) * rand(15, 10);
r5 = -20 + (-15 - (-20)) * rand(10, 10);
r6 = -0.7 + (0.1 - (-0.7)) * rand(50, 50);

% c
t1b = r2(4);

% BONUS
t1b = max(r2);

% d i
t1dd = r2([2, 3]);

% d ii
% The command r2(a, b) requires me to enter the number of rows and column
% of a matrix named r2. Since, r2 isn't a matrix, that command will return
% an invalid dimension error.

% e
t1e = r3([10, 12, 15]);

% f
t1f = r4(2:9);

% g
t1g = r4([1, 10:end]);

% h
% I can refer to the second-to-last element using end-1 
% e.g. r4(end-1)

% i
t1i = [r1(end), r2(end-1), r3(end-2)]

%% Exercise 5
% Task 2
clearvars
close all
clc

rng(1)
Y = randi(9, [7, 5]);

testi1 = Y(2, 1);            % (2,1) entry value
testi2 = Y(1, 2);
testi3 = Y([1,2], 1);        % (1,1), (2,1) entry value
testi4 = Y(5, [3,4,5]);    
testi5 = Y([1,2], [1,3]);
testi6 = Y([1,2], 1:2:5);    % (1,1), (1,3), 1, 5), ..., (2, 5) entry value
testi7 = Y(:, 1);            % first col
testi8 = Y(end, :);          % last row

% a
t2a = Y(3, 4);

% b
t2b = Y(4, :);

% c
t2c = Y(:, end);

% d
t2d = Y([2,6], [3, 5]);

% e
t2e = Y([end-1, end], [end-1, end]);

%% Exercise 5
% Task 3
clearvars
close all
clc

% a
rng(1);
r1 = rand(1, 10);
% r1(15)

% "Index exceeds the number of array elements. Index must not exceed 10.",
% was the error that was displayed. It indicates that the index that I'm 
% trying to retrieve a value from doesn't exist, because there are no more
% than 10 indices.

% b
% The command r2(a, b) requires me to enter the number of rows and column
% of a matrix named r2. Since, r2 isn't a matrix, that command will return
% an invalid dimension error.

% c
% sin = 359
% a = sin(1)
% b = sin(2)

% Even though sin is a mathematical function in MATLAB, we can use the word
% 'sin' as a variable name. Once we use it as a varible name and set
% sin=359, it becomes a 1-dimensional vector. That's why when sin(1)
% doesn't return any error instead it returns the value of index 1. On the
% other hand, sin(2) tries to get the value from index 2 and since there's
% no other indices than 1, it returns an IndexExeeded error.

%% Exercise 5
% Task 4
clearvars
close all
clc

% a
% --- i ---
A1 = [1, 1, -4; 2, -1, 2; 2, 1, 2];
b1 = [-2; 0; -6];
rref1 = rref([A1, b1])

% --- ii ---
A2 = [0, 1, -4; 2, 0, 2; 2, 2, -6];
b2 = [-3; 0; -6];
rref2 = rref([A2, b2])

% --- iii ---
A3 = [0, 1, -4; 2, 0, 2; 2, 2, -6];
b3 = [-3; 0; -5];
rref3 = rref([A3, b3])

% --- iv ---
A4 = [3, 0, -1, 4; 0, 2, -2, -1; 8, 0, 0, -2];
b4 = [0; 0; 0];
rref4 = rref([A4, b4])

% b
% --- i ---
% rref1 =
%      1.0000         0         0   -1.0000
%           0    1.0000         0   -3.0000
%           0         0    1.0000   -0.5000
% This matrix shows that there is a pivot in every column corresponding to
% a variable (x, y, z). Therefore, this system has a unique solution

% --- ii ---
% rref2 =
%      1     0     1     0
%      0     1    -4    -3
%      0     0     0     0
% The third row is all zeros, indicating there are infinitely many 
% solutions. The system has a free variable x3.

% --- iii ---
% rref3 =
% 
%      1     0     1     0
%      0     1    -4     0
%      0     0     0     1
% The last row, where 0 = 1, is a contradiction. Therefore, this system 
% has no solution.

% --- iv ---
% rref4 =
% 
%     1.0000         0         0   -0.2500         0
%          0    1.0000         0   -5.2500         0
%          0         0    1.0000   -4.7500         0
% There are 4 variable column but 3 rows, which means there are not enough
% equations to solve the system of equations. So, there are infinitely 
% many solutions.

% c
% Therefore, system i has unique solutions. Here is the solution vector:

% --- i ---
x1 = [rref1(:, end)];
A1*x1 == b1

%% Exercise 5
% Task 5
close all

% a
det_A1 = det(A1)
% Since the determinant is non-zero, the system has a unique solution.

det_A2 = det(A2)
% The determinant is zero, indicating that the system may have infinitely 
% many solutions or no solution.

det_A3 = det(A3)
% The determinant is zero, indicating that the system may have infinitely 
% many solutions or no solution.

% det_A4 = det(A4)
% Cannot obtain the determinant since the matrix isn't square in size.

% b
inv_A1 = inv(A1)

% inv_A2 = inv(A2)
% inv_A3 = inv(A3)
% inv_A4 = inv(A4)
% Since the determinant of A2 and A3 is zero, their inverse cannot be
% determined. It will return an infinity matrix and say that "Warning:
% Matrix is singular to working precision." because a zero division occurs
% while calculating the inverse of those matrices.
% A4 is not a square matrix, so inverse cannot be calculated either.

%% Exercise 5
% Task 6
clearvars
close all
clc

rng(3);
r = 2 * rand(1, 10)
ind = find(r > 1.2)
t6d = r(ind)
len_t6d = length(t6d)

%% Exercise 5
% Task 7
clearvars
close all
clc

% a
A0 = 7;
lambda = 0.4;
t = 0:0.1:10;

A_t = A0 * exp(-lambda * t);

plot(t, A_t);

% b
amount_remaining = A_t(end)

% c
% index is 51 at t=5
% The correponding amount of substance left is 0.947346982656289

% d
% Closest value to 4 is 3.998463446941704 at index 15
% The time value at index 15 is 1.4

% e
lambda1 = 0.2; % Smaller lambda
lambda2 = 0.6; % Larger lambda

A_t1 = A0 * exp(-lambda1 * t);
A_t2 = A0 * exp(-lambda2 * t);

hold on;
plot(t, A_t1, 'r');
plot(t, A_t2, 'g');
hold off;

% f
t_index_5_6 = find(t == 5.6);
t7f = [A_t(t_index_5_6), A_t1(t_index_5_6), A_t2(t_index_5_6)]

% g
B_t = A0 - A_t;
hold on;
plot(t, B_t, 'm--');

% h
title('Amount of Substance A Over Time');
xlabel('Time (t)');
ylabel('Amount of Substance A (A(t))');
grid on;

% i
index_below_threshold = find(A_t < 0.5, 1);
time_below_threshold = t(index_below_threshold)

% j
plot(5.6, A_t(find(t == 5.6)), 'bo', 'MarkerFaceColor', 'b');
plot(5.6, A_t1(find(t == 5.6)), 'ro', 'MarkerFaceColor', 'r');
plot(5.6, A_t2(find(t == 5.6)), 'go', 'MarkerFaceColor', 'g');
plot(5.6, B_t(find(t == 5.6)), 'mo', 'MarkerFaceColor', 'm');
legend('lambda = 0.4', 'lambda = 0.2', 'lambda = 0.6', 'The amount of B', ...
    'at t = 5.6', 'at t = 5.6', 'at t = 5.6', 'at t = 5.6');

