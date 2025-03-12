% Define the function
function [addition, subtraction, multiplication, division] = ex3p1(input_array)
    % Define the generated array [1, 2, 3, 4, 5]
    generated_array = [1, 2, 3, 4, 5];
    
    % Perform element-wise operations
    addition = input_array + generated_array;
    subtraction = input_array - generated_array;
    multiplication = input_array .* generated_array;
    division = input_array ./ generated_array;
    
    % Plot the results
    figure;
    hold on;
    
    % Plot addition
    plot(addition, '-o', 'DisplayName', 'Addition');
    
    % Plot subtraction
    plot(subtraction, '-x', 'DisplayName', 'Subtraction');
    
    % Plot multiplication
    plot(multiplication, '-s', 'DisplayName', 'Multiplication');
    
    % Plot division
    plot(division, '-d', 'DisplayName', 'Division');
    
    % Add title and labels
    title('Element-wise Operations');
    xlabel('Index');
    ylabel('Value');
    
    % Add legend
    legend show;
    grid on;
    hold off;
end

% Call the function with the input array
input_array = [10, 20, 30, 40, 50];
[add, sub, mul, div] = ex3p1(input_array);