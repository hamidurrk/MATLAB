function [X, Y, Z] = ex4p5(x, y)
    x = linspace(x(1), x(2), 10);
    y = linspace(y(1), y(2), 10);
    
    [X, Y] = meshgrid(x, y);
    
    Z = threehump(X, Y);
    
    figure;
    surf(X, Y, Z);
    title('Three-Hump Camel Function');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    
    set(gca, 'ZScale', 'log');
    
end

function z = threehump(x, y)    
    z = 2*x.^2 - 1.05*x.^4 + (x.^6)/6 + x.*y + y.^2;
end
