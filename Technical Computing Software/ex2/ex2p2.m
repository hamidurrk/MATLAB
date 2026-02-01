function result = ex2p2(height, mass)
    bmi = mass / (height^2);

    if bmi < 18.5
        result = "underweight";
    elseif bmi < 25
        result = "normal weight";
    else
        result = "overweight";
    end
end