% Exercise 6
% Task 1

% a
n = 10;
v = []; 
for s = 1:n
    v(s) = s^2;
end
v

n = 10;
v = []; 
for s = 1:n
    v(end + 1) = s^2;
end
v

n = 10;
v = zeros(1, n); 
for s = 1:n
    v(s) = s^2;
end
v

% b
n = 10;
v = []; 
for s = 1:n
    v = [v; 2*s+1];
end
v

n = 10;
v = []; 
for s = 1:n
    v(end + 1) = 2*s+1;
end
v = transpose(v);
v

n = 10;
v = zeros(1, n); 
for s = 1:n
    v(s) = 2*s+1;
end
v = transpose(v);
v

% c
n = 10;
v = []; 
for s = 1:2:n
    v(s) = s^2;
end
v
% The values will be set to zero by default.

% d
n = 100000;

tic
v = []; 
for s = 1:n
    v(s) = s^2;
end
toc 

tic
v = []; 
for s = 1:n
    v(end + 1) = s^2;
end
toc

tic
v = zeros(1, n); 
for s = 1:n
    v(s) = s^2;
end
toc

% e
n = 100000;
tic
v = []; 
for s = 1:n
    v = [v, 2*s];
end
toc

%% Exercise 6
% Task 2
clearvars
close all
clc

% a
n = 5;
if n == 0
    result = 1;
else
    result = 1;
    for i = 1:n
        result = result * i;
    end
end
disp(["Factorial of " + num2str(n)])
disp(result)

%% b
clearvars

n = 5;
if n == 0
    result = 1;
else
    result = 1;
    i = 1;
    while i <= n
        result = result * i;
        i = i + 1;
    end
end
disp(["Factorial of " + num2str(n)])
disp(result)

