% Exercise 4
% Task 1

% b
q = input('Enter the common ratio q: ');
n = input('Enter the number of terms n: ');

sum = 0;

for i = 0:n
    sum = sum + q^i; 
end

fprintf('The geometric sum of first n+1 terms is: %.2f\n', sum);

%% Exercise 4
% Task 2
clearvars
close all
clc

% a
for i = 3:3:27
    fprintf("%d\n", i)
    pause(0.5)
end

% c
for i = 1:9
    fprintf("%d\n", i*3)
    pause(0.5)
end

%% Exercise 4
% Task 3
clearvars
close all
clc

% a
tic
a = ['m', 'b', 'g', 'y', 'r'];
hold on
axis equal
axis([-10, 10, 0, 10])
for i = 5:9
    x = -i:0.1:i;
    y = sqrt(i^2 - x.^2);
    plot(x, y, a(i - 4))
end
xlabel('x')
ylabel('y')
title('Rainbow')
legend({'$y = \sqrt{i^2 - x^2}$'}, 'interpreter', 'latex')
toc

% b
% It takes around 0.33 plus minus 0.02 seconds to run the program

% c
% After removing the unnecessary outputs, it takes around 0.26 seconds

% d
% Now it takes only about 0.09 seconds


%% Exercise 4
% Task 4
clearvars
close all
clc

% a, b, c, d
x = -pi:0.1:2*pi;
y = sin(2*x);
plot(x, y, 'k');
title('Sine Curve')
hold on
% green star
for i = 1:length(x)
    if y(i) > 0.6
        plot(x(i), y(i), 'r*')
    elseif y(i) < -0.6
        plot(x(i), y(i), 'b*')
    else
        plot(x(i), y(i), 'g*')
    end
end


%% Exercise 4
% Task 5
clearvars
close all
clc

% a, b
v = [1, 2, 1, -1, 0, 5, -2, 3]; 
minimum = inf;
maximum = -inf;

for element = v
    if element < minimum 
        minimum = element; 
    end
    if element > maximum 
        maximum = element;
    end
end

fprintf("The smallest element in the vector is %d\n", minimum)
fprintf("The largest element in the vector is %d\n", maximum)
%% Exercise 4
% Task 6
clearvars
close all
clc

matrix = input('Enter the matrix (e.g., [1, 4, 5; 2, 3, 6]): ');
count_greater_than_3 = 0;
[rows, cols] = size(matrix);

for i = 1:rows
    for j = 1:cols
        if matrix(i, j) > 3
            count_greater_than_3 = count_greater_than_3 + 1; 
        end
    end
end

fprintf("The number of elements greater than 3 is: %d\n", count_greater_than_3)
