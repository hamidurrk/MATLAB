function out = ex2p5(op, a, b)
    switch op
        case '+'
            out = a + b;
        case '-'
            out = a - b;
        case '*'
            out = a .* b;
        case '/'
            out = a ./ b;
        otherwise
            out = NaN;
    end
end