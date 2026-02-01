%% MATLAB 1

function plot_on_ball(f,a,r)
    n = 300;
    rho = linspace(0,r,n);
    theta = linspace(0,2*pi,n);
    [R,T] = meshgrid(rho,theta);
    
    X = a(1) + R.*cos(T);
    Y = a(2) + R.*sin(T);
    Z = f(X,Y);
    
    surf(X,Y,Z,'EdgeColor','none')
    xlabel('x'), ylabel('y'), zlabel('f(x,y)')
    view(45,30)
end

f = @(x,y) x.*y;
a = [1 1];
r = 5;

plot_on_ball(f,a,r)

%% MATLAB 2

theta = linspace(0, 2*pi, 1000);
p_values = [linspace(1, 10, 50), linspace(11, 100, 20)];

figure;
h = plot(0, 0, 'LineWidth', 2);
axis equal;
grid on;
xlabel('x_1'); ylabel('x_2');

for p = p_values
    % p-norm: (|x1|^p + |x2|^p)^(1/p) = 1
    % using parametric form: x = cos(t) / (|cos(t)|^p + |sin(t)|^p)^(1/p)
    
    denom = (abs(cos(theta)).^p + abs(sin(theta)).^p).^(1/p);
    x1 = cos(theta) ./ denom;
    x2 = sin(theta) ./ denom;
    
    set(h, 'XData', x1, 'YData', x2);
    title(['p-norm unit circle, p = ', num2str(p, '%.2f')]);
    axis([-1.5 1.5 -1.5 1.5]);
    drawnow;
end