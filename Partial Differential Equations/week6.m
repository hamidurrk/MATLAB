%% MATLAB 1
clc; clear;

f = @(x,y) x.^2 - y.^2;

a = [0,0];
r = 1;

theta = linspace(0,2*pi,100000);
x = a(1) + r*cos(theta);
y = a(2) + r*sin(theta);

u_boundary = f(x,y);
mean_boundary = trapz(theta, u_boundary) / (2*pi);

u_center = f(a(1), a(2));

disp(mean_boundary)
disp(u_center)

%% MATLAB 2
clc; clear;

f = @(x) abs(x);

x = linspace(-1,1,2000);
eps_vals = [0.2 0.1 0.05];

figure
hold on

plot(x, f(x), 'k', 'LineWidth', 2)

for e = eps_vals
    eta = @(t) (1/e) * exp(-1./(1 - (t/e).^2)) .* (abs(t) < e);
    
    fe = zeros(size(x));
    for i = 1:length(x)
        y = linspace(-1,1,2000);
        fe(i) = trapz(y, f(y) .* eta(x(i) - y));
    end
    
    plot(x, fe, 'LineWidth', 1.5)
end

legend('f(x)=|x|', '\epsilon=0.2', '\epsilon=0.1', '\epsilon=0.05')
hold off
