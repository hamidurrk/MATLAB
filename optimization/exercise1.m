%% Problem 5
clc; clear; close all;

cities = [  0   0;      % x1 home
            0   1;      % x2
            3   5;      % x3
            -3  2;      % x4
            1   -4;     % x5
            -3  -3];    % x6

n = size(cities,1);

P = perms(2:n);

bestLength = inf;
bestRoute = [];

for i = 1:size(P,1)
    route = [1 P(i,:) 1]; % start and end at home
    L = 0;
    for j = 1:length(route)-1
        a = cities(route(j),:);
        b = cities(route(j+1),:);
        L = L + norm(a-b);
    end
    if L < bestLength
        bestLength = L;
        bestRoute = route;
    end
end

fprintf('Shortest route: ');
fprintf('%d ', bestRoute);
fprintf('\nShortest length = %.4f\n', bestLength);

figure;
hold on; grid on; axis equal;
plot(cities(:,1), cities(:,2), 'ko', 'MarkerFaceColor','y', 'MarkerSize',8);
plot(cities(bestRoute,1), cities(bestRoute,2), 'r-o', 'LineWidth',2);
text(cities(:,1)+0.1, cities(:,2)+0.1, {'x1','x2','x3','x4','x5','x6'});
title('Shortest Traveling Salesperson Route');
xlabel('x'); ylabel('y');

%% Problem 6
clc; clear; close all;

f   = @(x) x.*sin(2*x);
fp  = @(x) sin(2*x) + 2*x.*cos(2*x);        % f'(x)
fpp = @(x) 4*cos(2*x) - 4*x.*sin(2*x);      % f''(x)

starts = linspace(-pi,pi,30);
roots = [];

for s = starts
    x = s;
    for k = 1:40              
        x = x - fp(x)/fpp(x);
    end
    if abs(fp(x)) < 1e-6 && x>=-pi && x<=pi
        roots(end+1) = x; 
    end
end

roots = unique(round(roots,6));

candidates = [roots, -pi, pi];
values = f(candidates);

[fmin,idx] = min(values);
xmin = candidates(idx);

fprintf('Stationary points: '); fprintf('%.6f  ', roots); fprintf('\n');
fprintf('Global minimizer: x = %.6f\n', xmin);
fprintf('Minimum value: f(x) = %.6f\n', fmin);

xx = linspace(-pi,pi,2000);
figure; grid on; hold on;
plot(xx,f(xx),'LineWidth',2);
plot(xmin,fmin,'ro','MarkerFaceColor','r');
title('Minimization of f(x) = x sin(2x) using Newton''s Method');
xlabel('x'); ylabel('f(x)');
