function tf = ex2p4(year)
    if (mod(year,4) == 0) && ( (mod(year,100) ~= 0) || (mod(year,400) == 0) )
        tf = true;
    else
        tf = false;
    end
end