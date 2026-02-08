function mse = ex4p3()
    data = load('temperature_data.mat');
    
    x = data.date(:);  % Ensure column vector
    y = data.temperature(:);  % Ensure column vector
    
    X = [ones(size(x)), x, x.^2];
    
    b = X \ y;
    
    y_fitted = X * b;
    
    error_vector = y - y_fitted;
    mse = mean(error_vector.^2);
    
    figure;
    plot(x, y, '*'); % Original data with '*'
    hold on;
    plot(x, y_fitted, '--', 'LineWidth', 2); % Fitted curve with '--'
    hold off;
    
    xlabel('Date');
    ylabel('Temperature');
    title(sprintf('fitted model, MSE = %.4f', mse));
    legend('Original Data', 'Fitted Curve', 'Location', 'best');
    grid on;
end
