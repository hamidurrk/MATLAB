%% Task 4
clc; clear;

f = @(x) cos(x);
I = sin(1);

ns = [20 40 60 80 100];
err = zeros(size(ns));
sn  = zeros(size(ns));

for k = 1:numel(ns)
    n = ns(k);
    x = (0:n)/n;      
    y = f(x);
    h = 1/n;

    sn(k) = sum( h * (y(2:end) + y(1:end-1)) / 2 );

    err(k) = abs(I - sn(k));
end

T = table(ns(:), sn(:), err(:), 'VariableNames', {'n','s_n','abs_error'});
disp(T)

%% Task 5
clc; clear;

f = @(x) cos(x);
I = sin(1);

ns = [20 40 60 80 100];          
err = zeros(size(ns));
sn  = zeros(size(ns));

for k = 1:numel(ns)
    n = ns(k);
    h = 1/n;
    x = 0:h:1;                   
    y = f(x);

    sn(k) = (h/3) * ( y(1) + y(end) ...
                    + 4*sum(y(2:2:end-1)) ...
                    + 2*sum(y(3:2:end-2)) );

    err(k) = abs(I - sn(k));
end

T = table(ns(:), sn(:), err(:), 'VariableNames', {'n','s_n','abs_error'});
disp(T)