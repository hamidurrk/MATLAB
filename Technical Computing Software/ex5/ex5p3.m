function table = ex5p3(rows, cols)
    table = zeros(rows, cols); 

    for i = 1:rows
        for j = 1:cols
            table(i, j) = i * j;
        end
    end

end