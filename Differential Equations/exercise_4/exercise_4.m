%% Exercise 4
% Task 1
clearvars
clc
close all

% b
odefun = @(t,y) [y(2); -y(1)];
tspan = [0 10];
y0 = [1; 0];
[T, Y] = ode45(odefun, tspan, y0);

figure;
plot(T, Y(:,1), 'b-');
hold on;

syms y(t)
Dy = diff(y);
ode = diff(y,t,2) + y == 0;
cond = [y(0)==1, Dy(0)==0];
ySol(t) = dsolve(ode, cond);

ySym = double(ySol(T));
plot(T, ySym, 'ro');

title("1(b) Solution of y'' + y = 0");
xlabel('t');
ylabel('y(t)');
legend('Numerical', 'Symbolic (dsolve)');

disp('Sizes of T and Y:');
disp(['Size of T: ', num2str(size(T))]);
disp(['Size of Y: ', num2str(size(Y))]);


%% c
clearvars
clc
close all

% 1(c)
odefun_c = @(t,y) [y(2); -y(1) - 2*y(2)];
y0_c = [1; 0];
[t_c, y_c] = ode45(odefun_c, [0 10], y0_c);

figure;
plot(t_c, y_c(:,1));
title('1(c) Solution of y"" + 2y'' + y = 0');
xlabel('t');
ylabel('y(t)');

%% Exercise 4
% Task 2
clearvars
clc
close all

% a
syms y(t);
eqn = diff(y,t,2) + y == 0;
V = odeToVectorField(eqn);
odefun_symbolic = matlabFunction(V, 'Vars', {'t','Y'});

disp('Generated function for 2(a):');
disp(odefun_symbolic);

% b

[T_sym, Y_sym] = ode45(odefun_symbolic, [0 10], [1; 0]);

figure;
plot(T_sym, Y_sym(:,1), 'b-');
title('2(b) Solution using Symbolic Conversion');
xlabel('t');
ylabel('y(t)');
legend('Numerical');


%% c
clearvars
clc 
close all

syms y(t);
eqn = diff(y,t,3) + 2*diff(y,t) + y == 0;
V = odeToVectorField(eqn);
odefun_3rd = matlabFunction(V, 'Vars', {'t','Y'});

y0_3rd = [1; 1; 0];
[T_3rd, Y_3rd] = ode45(odefun_3rd, [0 10], y0_3rd);

figure;
plot(T_3rd, Y_3rd(:,1));
title('2(c) Solution of y^{(3)} + 2y'' + y = 0');
xlabel('t');
ylabel('y(t)');


%% d i
clearvars
clc
close all

syms y(t);
eqn_i = t*diff(y,t,2) + t^2*diff(y,t) + y == 0;
try
    V_i = odeToVectorField(eqn_i);
    odefun_i = matlabFunction(V_i, 'Vars', {'t','Y'});
    [T_i, Y_i] = ode45(odefun_i, [1 5], [1; 0]);
    
    figure;
    plot(T_i, Y_i(:,1), 'b-');
    title('2(d)(i) Numerical Solution of t y"" + t^2 y'' + y = 0');
    xlabel('t');
    ylabel('y(t)');
catch
    disp('Symbolic conversion failed for 2(d)(i)');
end

% d ii
syms y(t);
eqn_ii = t^2*diff(y,t,2) + t*diff(y,t) + y == 0;
V_ii = odeToVectorField(eqn_ii);
odefun_ii = matlabFunction(V_ii, 'Vars', {'t','Y'});

% Numerical solution
[T_ii, Y_ii] = ode45(odefun_ii, [1 10], [1; 0]);

% Symbolic solution
dy = diff(y);
cond_ii = [y(1) == 1, dy(1) == 0];
ySol_ii(t) = dsolve(eqn_ii, cond_ii);
ySym_ii = double(ySol_ii(T_ii));

figure;
plot(T_ii, Y_ii(:,1), 'b-');
hold on;
plot(T_ii, ySym_ii, 'ro');
title('2(d)(ii) Solution of t^2 y"" + t y'' + y = 0');
xlabel('t');
ylabel('y(t)');
legend('Numerical', 'Symbolic');

%% Exercise 4
% Task 3
clearvars
clc
close all

m = 1;
k = 0.5;
g = 9.81;

tspan = [0, 10];  

y0 = [100; 0];

% Without air resistance
f_noDrag = @(t,y) [ y(2); -g ];           

% With air resistance
f_withDrag = @(t,y) [ y(2); -g - (k/m)*y(2) ];

[t_noDrag, y_noDrag]     = ode45(f_noDrag,    tspan, y0);
[t_withDrag, y_withDrag] = ode45(f_withDrag,  tspan, y0);

% position vs. time plot
figure('Name','Falling Object - Position vs Time','NumberTitle','off');
hold on; grid on;
plot(t_noDrag,    y_noDrag(:,1),    'b-','LineWidth',2, 'DisplayName','No Drag');
plot(t_withDrag,  y_withDrag(:,1),  'r-','LineWidth',2, 'DisplayName','With Drag');
xlabel('Time (s)');
ylabel('Position s(t) (m)');
title('Falling Object with/without Air Resistance');
legend('Location','Best');
hold off;

%% Animation

figure('Name','Falling Object Animation','NumberTitle','off');
hold on; grid on;
axis([tspan(1), tspan(2), 0, 110]);
xlabel('Time (s)');
ylabel('Height (m)');
title('Animation: Falling Object (Upward = Positive)');

plot(t_noDrag,    y_noDrag(:,1),    'b--','HandleVisibility','off');
plot(t_withDrag,  y_withDrag(:,1),  'r--','HandleVisibility','off');

markerNoDrag   = plot(t_noDrag(1),   y_noDrag(1,1),   'bo','MarkerFaceColor','b','MarkerSize',8,...
                      'DisplayName','No Drag');
markerWithDrag = plot(t_withDrag(1), y_withDrag(1,1), 'ro','MarkerFaceColor','r','MarkerSize',8,...
                      'DisplayName','With Drag');
legend('No Drag','With Drag','Location','Best');

idxNoDrag   = 1;  
idxWithDrag = 1;  

tCurrent = 0;     
tEnd     = max(t_noDrag(end), t_withDrag(end));

tStep = 0.01;  

while tCurrent <= tEnd
    sNoDrag   = interp1(t_noDrag,   y_noDrag(:,1),   tCurrent, 'linear','extrap');
    sWithDrag = interp1(t_withDrag, y_withDrag(:,1), tCurrent, 'linear','extrap');
    
    set(markerNoDrag,   'XData', tCurrent, 'YData', sNoDrag);
    set(markerWithDrag, 'XData', tCurrent, 'YData', sWithDrag);
    
    pause(tStep);
    
    tCurrent = tCurrent + tStep;
end

hold off;
