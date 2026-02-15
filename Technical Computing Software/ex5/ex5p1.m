function result = ex5p1(vec)
    result = false(size(vec));

    for i = 1:numel(vec)
        result(i) = prime(vec(i));
    end

end