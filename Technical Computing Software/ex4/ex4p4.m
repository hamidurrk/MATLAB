function distance = ex4p4()
    data_struct = load('h03t2data.mat');
    signal = data_struct.data;
    
    fs = 1000; % Hz (1 kHz)
    
    control = signal(1, :);
    x = signal(2, :);
    y = signal(3, :);
    
    threshold = 1.0; % Value above which device is considered "on"
    is_on = control > threshold;
    
    on_indices = find(is_on);
    
    if isempty(on_indices)
        error('No powered-on state detected');
    end
    
    start_idx = on_indices(1);
    end_idx = on_indices(end);
    
    t = (on_indices - start_idx) / fs; % Time in seconds
    x_on = x(on_indices);
    y_on = y(on_indices);
    
    x_start = x_on(1);
    y_start = y_on(1);
    x_end = x_on(end);
    y_end = y_on(end);
    
    distance = sqrt((x_end - x_start)^2 + (y_end - y_start)^2);
    
    figure;
    
    subplot(1, 2, 1);
    plot3(t, x_on, y_on, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel('X Position (mm)');
    zlabel('Y Position (mm)');
    title(sprintf('Distance traveled %.2f mm', distance));
    view(45, 30); % Set a good viewing angle
    
    hold on;
    plot3(t(1), x_start, y_start, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot3(t(end), x_end, y_end, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    legend('Path', 'Start', 'End', 'Location', 'best');
    hold off;
    
    subplot(1, 2, 2);
    plot(x_on, y_on, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('X Position (mm)');
    ylabel('Y Position (mm)');
    title('X-Y Path');
    axis equal;
    
    hold on;
    plot(x_start, y_start, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot(x_end, y_end, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    legend('Path', 'Start', 'End', 'Location', 'best');
    hold off;
end
