clearvars; close all;

% ---- Problem data ----
A = 3; B = -2; C = 3;
f  = @(x,y) A*x.^2 - B*x.*y + C*y.^2 + x - y;  % objective function

% forward-difference approximations
h = 1e-5;
dfdx_fd = @(x,y) (f(x+h,y) - f(x,y)) / h;
dfdy_fd = @(x,y) (f(x,y+h) - f(x,y)) / h;

% ---- Newton with dynamic alpha ----
max_iter  = 100;
tol       = 1e-5;
x = 5; y = 7;         % starting point

history_xy = [x; y];
history_f  = f(x,y);

for k = 1:max_iter
    %--- 1) build numerical gradient and normalize
    gx = dfdx_fd(x,y);
    gy = dfdy_fd(x,y);
    g = [gx; gy];
    d = g / norm(g);
    
    %--- 2) form phi'(alpha)
    phi_prime = @(alpha) ...
      dfdx_fd( x - alpha*d(1), y - alpha*d(2) ) * d(1) + ...
      dfdy_fd( x - alpha*d(1), y - alpha*d(2) ) * d(2);
    
    %--- 3) bracket until sign change for phi'
    a_lo = 0;  f_lo = phi_prime(a_lo);
    a_hi = 1;  f_hi = phi_prime(a_hi);
    while f_lo * f_hi > 0
        a_hi = 2*a_hi;
        f_hi = phi_prime(a_hi);
        if a_hi > 1e6
            error('Could not bracket phi'' root.');
        end
    end

    %--- 4) solve phi'(alpha) = 0
    alpha_opt = fzero(phi_prime, [a_lo, a_hi]);
    
    %--- 5) update (x,y)
    x_new = x - alpha_opt * d(1);
    y_new = y - alpha_opt * d(2);
    
    % record history
    history_xy(:,end+1) = [x_new; y_new];
    history_f(end+1)    = f(x_new,y_new);
    
    % check convergence
    if norm([x_new - x; y_new - y]) < tol
        fprintf('Converged after %d iterations.\n', k);
        break
    end
    
    % if norm([x_new + 0.25; y_new - 0.25]) < tol
    %     fprintf('Reached the desired point (-0.25, 0.25) after %d iterations.\n', k);
    %     break
    % end
    
    % prepare for next iteration
    x = x_new;
    y = y_new;
end

% ---- Visualization ----
% surface + path
[xg,yg] = meshgrid(-2:.1:8,-2:.1:8);
zg = f(xg,yg);
figure;
surf(xg,yg,zg, 'DisplayName', 'Function Surface'); hold on;

% Plot initial and final points
plot3(history_xy(1,1), history_xy(2,1), history_f(1), 'ro', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Initial Point');
plot3(history_xy(1,end), history_xy(2,end), history_f(end), 'go', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Final Point');

xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('Newton path with numerical derivatives and dynamic \alpha');
view(-30, 30); 
colorbar;
grid on;
legend('Location', 'best');

% contour + path
figure;
contour(xg,yg,zg,50); hold on;
% Plot path line
plot(history_xy(1,:), history_xy(2,:), 'b-','LineWidth',2);

% Plot initial and final points
plot(history_xy(1,1), history_xy(2,1), 'ro', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Initial Point');
plot(history_xy(1,end), history_xy(2,end), 'go', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Final Point');

xlabel('x'); ylabel('y');
title('Contour + Newton path');
grid on;
legend('Location', 'northwest');

% Create table data
n_rows = length(history_xy(1,:));
table_data = cell(n_rows, 7);

% Fill initial row (k=0)
table_data{1,1} = 0;
table_data{1,2} = sprintf('(%.5f, %.5f)', history_xy(1,1), history_xy(2,1));
g0 = [dfdx_fd(history_xy(1,1), history_xy(2,1)); 
      dfdy_fd(history_xy(1,1), history_xy(2,1))];
table_data{1,3} = sprintf('(%.5f, %.5f)', g0(1), g0(2));
d0 = g0 / norm(g0);
table_data{1,4} = sprintf('(%.5f, %.5f)', d0(1), d0(2));
table_data{1,5} = '-';
table_data{1,6} = '-';
table_data{1,7} = '-';

% Fill remaining rows
for i = 2:n_rows
    table_data{i,1} = i-1;
    table_data{i,2} = sprintf('(%.5f, %.5f)', history_xy(1,i-1), history_xy(2,i-1));
    
    % Calculate gradient at current point
    g = [dfdx_fd(history_xy(1,i-1), history_xy(2,i-1)); 
         dfdy_fd(history_xy(1,i-1), history_xy(2,i-1))];
    table_data{i,3} = sprintf('(%.5f, %.5f)', g(1), g(2));
    
    % Normalized gradient
    d = g / norm(g);
    table_data{i,4} = sprintf('(%.5f, %.5f)', d(1), d(2));
    
    % For the last row, calculate what the next point would be
    if i == n_rows
        % Form phi'(alpha)
        phi_prime = @(alpha) ...
            dfdx_fd(history_xy(1,i-1) - alpha*d(1), history_xy(2,i-1) - alpha*d(2)) * d(1) + ...
            dfdy_fd(history_xy(1,i-1) - alpha*d(1), history_xy(2,i-1) - alpha*d(2)) * d(2);
        
        % Bracket and solve for alpha
        a_lo = 0;  f_lo = phi_prime(a_lo);
        a_hi = 1;  f_hi = phi_prime(a_hi);
        while f_lo * f_hi > 0
            a_hi = 2*a_hi;
            f_hi = phi_prime(a_hi);
            if a_hi > 1e6
                break;
            end
        end
        alpha_opt = fzero(phi_prime, [a_lo, a_hi]);
        
        % Calculate theoretical next point
        x_next = history_xy(1,i-1) - alpha_opt * d(1);
        y_next = history_xy(2,i-1) - alpha_opt * d(2);
        
        table_data{i,5} = sprintf('%.5f', alpha_opt);
        table_data{i,6} = sprintf('(%.5f, %.5f)', x_next, y_next);
        table_data{i,7} = sprintf('%.4e', norm([x_next; y_next] - history_xy(:,i-1)));
    else
        diff = history_xy(:,i) - history_xy(:,i-1);
        alpha = norm(diff) / norm(d);
        table_data{i,5} = sprintf('%.5f', alpha);
        table_data{i,6} = sprintf('(%.5f, %.5f)', history_xy(1,i), history_xy(2,i));
        table_data{i,7} = sprintf('%.4e', norm(diff));
    end
end

% Display the table
fprintf('\nOptimization Progress:\n');
headers = {'Iteration', '(x,y)', '∇f(x,y)', '∇f(x,y) normalized', 'α', 'next(x,y)', 'Error'};
disp(array2table(table_data, 'VariableNames', headers));
