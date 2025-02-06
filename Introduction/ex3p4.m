function [value1, value2, value3] = ex3p4(A, B, C)
    value1 = vertcat([A, B], C);
    
    value2 = value1(2:4, :); 
    
    value3 = value1((value1 > 5) & (value1 < 15));
    value3 = value3(:); 
end


mat1 = randi(9, 3, 3)
mat2 = randi(9, 3, 3)
mat3 = randi(9, 3, 6)

value1 = vertcat([mat1, mat2], mat3)

value2 = value1(2:4, :)

value3 = value1((value1 > 5) & (value1 < 15))
value3 = value3(:)