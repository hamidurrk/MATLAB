function f = ex5p4(n)
    if n == 0
        f = 0;
    elseif n == 1
        f = 1;
    else
        f = ex5p4(n-1) + ex5p4(n-2);
    end

end