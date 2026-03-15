v0 = 45;
theta_deg = 20;
theta = deg2rad(theta_deg);

y0 = 2;
m = 0.145;
Cd = 0.3;
rho = 1.204;
d = 0.075;
g = 9.81;

A = pi*(d/2)^2;
k = 0.5*rho*Cd*A;

vx0 = v0*cos(theta);
vy0 = v0*sin(theta);

disp(['A = ', num2str(A)])
disp(['k = ', num2str(k)])
disp(['vx0 = ', num2str(vx0)])
disp(['vy0 = ', num2str(vy0)])