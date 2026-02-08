function [value1, value2, value3] = ex3p4(A, B, C)
    value1 = vertcat([A, B], C);
    
    value2 = value1(2:4, :); 
    
    value3 = value1((value1 > 5) & (value1 < 15));
    value3 = value3(:); 
end