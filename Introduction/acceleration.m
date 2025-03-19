file = load("./Introduction/acceleration.mat");
% stop_time = max(data.ts.Time)

plot(file.ts.Time, file.ts.Data)