function [elem_val, second_col, second_col_trans, product] = ex3p3(input_matrix)
    internal_matrix = [
        17, 24, 1, 8, 15;
        23, 5, 7, 14, 16;
        4, 6, 13, 20, 22;
        10, 12, 19, 21, 3;
        11, 18, 25, 2, 9
    ];
    
    elem_val = input_matrix(2, 4);
    
    second_col = input_matrix(:, 2);
    
    second_col_trans = second_col';
    
    product = input_matrix * internal_matrix;
end