%% MATLAB 1
clc; clear; close all;

f1 = @(v) abs(v(1)*v(2) - 4) + (v(1) - 7)^2 + (v(2) - 3)^2;
x0_1 = [0; 0];

opts_search = optimset('Display','off','TolX',1e-10,'TolFun',1e-10,'MaxIter',2e4,'MaxFunEvals',2e4);
[x_fminsearch, fval_fminsearch] = fminsearch(f1, x0_1, opts_search);

opts_unc = optimoptions('fminunc','Algorithm','quasi-newton','Display','off', ...
	'OptimalityTolerance',1e-10,'StepTolerance',1e-10,'FunctionTolerance',1e-10, ...
	'MaxIterations',2e4,'MaxFunctionEvaluations',2e4);
[x_fminunc, fval_fminunc] = fminunc(f1, x0_1, opts_unc);

fprintf('fminsearch: x = [%.10f, %.10f], f = %.10f\n', x_fminsearch(1), x_fminsearch(2), fval_fminsearch);
fprintf('fminunc:    x = [%.10f, %.10f], f = %.10f\n\n', x_fminunc(1), x_fminunc(2), fval_fminunc);

%% MATLAB 2
clc; clear; close all;

f2 = @(x) x(1)^2 + x(2)^2;
g = @(x) 1 - x(1) - x(2);
x0_2 = [0.2; 0.2];

opts2 = optimset('Display','off','TolX',1e-10,'TolFun',1e-10,'MaxIter',2e4,'MaxFunEvals',2e4);

r = 1;
xp = x0_2;
for k = 1:10
	phi_pen = @(x) f2(x) + r*max(0, g(x))^2;
	xp = fminsearch(phi_pen, xp, opts2);
	r = 10*r;
end
fprintf('2a (penalty): x = [%.10f, %.10f], f = %.10f, g = %.10f\n', xp(1), xp(2), f2(xp), g(xp));

mu = 1;
xb = [0.8; 0.8];
for k = 1:10
	phi_bar = @(x) f2(x) - mu*log(-g(x));
	xb = fminsearch(phi_bar, xb, opts2);
	if g(xb) >= 0
		xb = xb + 1e-3*[1; 1];
	end
	mu = 0.1*mu;
end
fprintf('2b (barrier): x = [%.10f, %.10f], f = %.10f, g = %.10f\n', xb(1), xb(2), f2(xb), g(xb));
