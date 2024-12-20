% Exercise 1
% Task 1
clearvars
close all
clc

sini = sin(pi/2)
kosini = cos(3*pi)

% Calculating sin(pi)

answer = sin(pi)

% Here the sin function is being called and the parameter is pi. This
% function is calulating the value of sin of pi radian or 180 degrees.
% The answer should be 0 but it's showing 1.2245e-16, which is not equal to
% 0 but a really small number. I think it is because pi here is the
% floating point approximation of the actual pi, (i.e. 3.1416) Also floating point data
% has some tolerances to consider and within the tolerance, this value can
% be considered 0. 

%% Exercise 1 
% Task 2
clearvars
close all
clc

T2i     = (23.4*2^2 - 238)/(9^2 + 34.2)
T2ii    = (8/9)^3*32 - 3^6/(7^7 - 555)
T2iii   = 4*sqrt(2200) - 45
T2iv    = sqrt(4^5) + log(exp(-5))
T2v     = sin(pi) - abs(cos(pi)/2)
T2vi    = (23 + exp(-2))/(exp(4) + log(1023))
T2vii   = (tan((pi/6)*log(8)))/(sqrt(17) + 2)
T2viii  = nthroot(-8, 3) + nthroot(6103515625, 7)
T2ix    = cos(5*pi/6)*(sin(8*pi/7))^2

%% Exercise 1 
% Task 4
clearvars
close all
clc

% a
ects_math = 5

% b
hours_corresponding_one_credit = 1600/60

% c 
hours_expected_from_math_student = ects_math * hours_corresponding_one_credit

% d
hours_needed_per_week = hours_expected_from_math_student/7
% Here, the exam week has been excluded from the calculation

% e
ects_test = 3
hours_corresponding_one_credit = 1600/60
hours_expected_from_student = ects_test * hours_corresponding_one_credit
hours_needed_per_week = hours_expected_from_student/7
% Again, the exam week has been excluded from the calculation


%% Exercise 1
% Task 5
clearvars
close all
clc

% b 
syms x;
f(x) = x^2
subs(f, x, 2)

% c 
fplot(f(x))
% Here, I have plotted the function f(x) using fplot function. Among the
% several parameters of fplot, including the function I want to plot,
% there's another parameter where it takes the interval of the x value to
% plot on the graph in the form of [a, b]. If the interval is omitted then
% the function uses it's default interval to plot the function. The default
% interval is [-5, 5]

% d
eq = x^2 - 1
solutions = solve(eq, x)

figure
fplot(eq)
hold on
fplot(solutions)
grid on

%% Exercise 1
% Task 6
clearvars
close all
clc

% a
% 2 = 2
% verified that 2 = 2 doesn't work. '=' is an assignment operator in
% MATLAB

% b 
2 == 2      % first equation
2 == 3      % second equation
% The error message from (a) described that '=' cannot be used to test the
% validity of an equation. Instead, it told me to use '==' operator for
% testing. We know that the first is True and the second one is False.
% Here in MATLAB, '==' has been used to test the validity. The first
% equation returns a logical 1 and the second one returns a logical 0,
% which is correct.

% c
syms x;
eq = x^2 - 1
solutions = solve(eq, x)

% bonus
syms x real
% first inequality
ineq1 = x^2 > 1
solution1 = solve(ineq1, x)
disp('Solution of the first inequality:')
disp(solution1)

% second inequality
ineq2 = abs(x - 5) > abs(x) + 1
solution2 = solve(ineq2, x)
disp('Solution of the second inequality:')
disp(solution2)