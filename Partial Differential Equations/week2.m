%% MATLAB 1

clear; clc;

% Grid parameters
Nx = 80;              % number of x steps
Nt = 200;             % number of time steps
x = linspace(0,1,Nx);
t = linspace(0,1,Nt);

dx = x(2) - x(1);
dt = t(2) - t(1);

% CFL number
lambda = 2 * dt / dx;
fprintf('CFL number = %.3f\n', lambda);

% Solution array
u = zeros(Nx, Nt);

% Initial condition
u(:,1) = 0;

% Time marching
for n = 1:Nt-1

    % Boundary condition at x = 0
    u(1,n) = 0.5 * atan(2 * t(n));

    for i = 2:Nx
        source = 1 / (x(i)^2 + 1);
        u(i,n+1) = u(i,n) ...
                   - 2 * dt/dx * (u(i,n) - u(i-1,n)) ...
                   + dt * source;
    end
end

% Exact solution
[X,T] = meshgrid(x,t);
u_exact = 0.5 * atan(X') - 0.5 * atan(X' - 2*T');

% Plot numerical vs exact at final time
figure;
plot(x, u(:,end), 'r-o', 'LineWidth', 1.2); hold on;
plot(x, u_exact(:,end), 'k--', 'LineWidth', 1.5);
legend('Numerical (Backward x, Forward t)', 'Exact');
xlabel('x');
ylabel('u(x,t)');
title('Comparison at t = 1');
grid on;

% Case 1: CFL ≤ 1 (Stable)
% Solution closely matches the exact solution
% Smooth propagation along characteristics
% Errors remain bounded
% 
% Case 2: CFL > 1 (Unstable)
% Oscillations appear
% Solution blows up rapidly
% Completely deviates from exact solution

%% MATLAB 2

clear; clc;

% Grid
Nx = 80;
Nt = 200;
x = linspace(0,1,Nx);
t = linspace(0,1,Nt);

dx = x(2) - x(1);
dt = t(2) - t(1);

% CFL number
lambda = 2 * dt / dx;
fprintf('CFL number = %.3f\n', lambda);

% Solution array
u = zeros(Nx, Nt);

% Initial condition
u(:,1) = 1 ./ (x.^2 + 1);

% Time stepping
for n = 1:Nt-1

    % Boundary condition at x = 1
    u(end,n) = 1 / ((1 + 2*t(n))^2 + 1);

    for i = 2:Nx
        u(i,n+1) = u(i,n) ...
                 + 2 * dt/dx * (u(i,n) - u(i-1,n));
    end
end

% Exact solution
[X,T] = meshgrid(x,t);
u_exact = 1 ./ ((X' + 2*T').^2 + 1);

% Plot at final time
figure;
plot(x, u(:,end), 'r-o', 'LineWidth', 1.2); hold on;
plot(x, u_exact(:,end), 'k--', 'LineWidth', 1.5);
legend('Numerical (Backward x, Forward t)', 'Exact');
xlabel('x');
ylabel('u(x,t)');
title('Comparison at t = 1');
grid on;


% It doesn't work. It becomes unstable
% 
% Using backward differences for u_x and forward differences for u_t does not 
% produce a stable solution for this problem.
% Although the spatial discretization respects the upwind direction, 
% the forward time discretization leads to numerical instability.
% The solution exhibits oscillations and blow-up due to violation of the 
% stability condition for this transport equation.

