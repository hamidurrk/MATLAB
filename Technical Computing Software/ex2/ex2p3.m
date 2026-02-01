function group = ex2p3(age, hasID)
    if (age <= 18) && (hasID)
        group = "youth";
    elseif (age > 18) && (age < 35)
        group = "open";
    elseif (age >= 35) && (hasID)
        group = "masters";
    else
        group = "open";
    end
end