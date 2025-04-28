clc; clearvars; close all;

% ---- Problem data ----
A = 3; B = -2; C = 3;
f  = @(x,y) A*x.^2 - B*x.*y + C*y.^2 + x - y;  % objective function

% forward-difference approximations
h = 1e-6;
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
plot3(history_xy(1,:), history_xy(2,:), history_f, 'r-o','LineWidth',2);
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('Newton path with numerical derivatives and dynamic \alpha');
view(-30, 30); 
colorbar;
grid on;

% contour + path
figure;
contour(xg,yg,zg,50); hold on;
plot(history_xy(1,:), history_xy(2,:), 'r-o','LineWidth',2);
xlabel('x'); ylabel('y');
title('Contour + Newton path');
grid on;
