function [T, t] = ex4p2()
    S = load("lpr_temp.mat");
    data = S.data;   % 145 x 3 matrix

    hours = data(:,1);
    minutes = data(:,2);
    temp = data(:,3);

    idx = (hours > 12 & hours < 17) | ...
          (hours == 12) | ...
          (hours == 17 & minutes == 0);

    T = temp(idx);

    year  = 2020 * ones(size(T));
    month = 6 * ones(size(T));
    day   = 1 * ones(size(T));

    t = datetime(year, month, day, hours(idx), minutes(idx), zeros(size(T)));

    figure
    plot(t, T, "-o")
    title("Temperature on 1st of July")
    xlabel("Time")
    ylabel("Temperature ^{o}C")
    grid on
end