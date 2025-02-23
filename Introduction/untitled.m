function [addition, subtraction, multiplication, division] = ex3p1(input_array)
    % Define the generated array [1, 2, 3, 4, 5]
    generated_array = [1, 2, 3, 4, 5];
    
    % Perform element-wise operations
    addition = input_array + generated_array;
    subtraction = input_array - generated_array;
    multiplication = input_array .* generated_array;
    division = input_array ./ generated_array;
end
input_array = [10, 20, 30, 40, 50];
[add, sub, mul, div] = ex3p1(input_array)