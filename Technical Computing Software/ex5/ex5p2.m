function nextPrime = ex5p2(n)
    candidate = n + 1;
    while ~prime(candidate)
        candidate = candidate + 1;
    end
    nextPrime = candidate;

end