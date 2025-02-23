%% Exercise 5
% Task 1
clearvars
clc
close all

% b

f = @(t,Y) [Y(2);
            -2*Y(2) - 3*Y(1)];

Y0 = [-2; 5];

tSpan = [1 5];
[tSol, YSol] = ode45(f, tSpan, Y0);


figure;
plot(tSol, YSol(:,1), 'b-o', 'DisplayName','y(t)'); hold on;
plot(tSol, YSol(:,2), 'r-x', 'DisplayName','y''(t)');
grid on;
xlabel('t');
ylabel('Solution');
title('Numerical Solution of y'''' + 2y'' + 3y = 0');
legend('Location','best');

%% c
clc; clear; close all;

f = @(t,XY) [ 0.3*XY(1) - 0.05*XY(1)*XY(2);
             -0.5*XY(2) + 0.003*XY(1)*XY(2) ];

XY0 = [20; 500];

tSpan = [0 200];
[tSol, XYSol] = ode45(f, tSpan, XY0);

xSol = XYSol(:,1);
ySol = XYSol(:,2);

figure;
plot(tSol, xSol, 'b-', 'LineWidth',2, 'DisplayName','x(t) - Prey'); hold on;
plot(tSol, ySol, 'r-', 'LineWidth',2, 'DisplayName','y(t) - Predator');
grid on;
xlabel('Time');
ylabel('Population');
title('Predator-Prey System (ESE6)');
legend('Location','best');

%% f
syms x(t) y(t)
eqn1 = diff(x,t) == 0.3*x - 0.05*x*y;
eqn2 = diff(y,t) == -0.5*y + 0.003*x*y;

cond = [ x(0)==20, y(0)==500 ];

[xSolSym, ySolSym] = dsolve(eqn1, eqn2, cond);

xSolSym, ySolSym

%% e
clc; clear; close all;

syms y1(t) y2(t)

ode1 = diff(y1, t) == y1 - 2*y2;
ode2 = diff(y2, t) == 3*y1 - 4*y2;

cond1 = y1(1) == -1;
cond2 = y2(1) == 1;

sol = dsolve(ode1, ode2, cond1, cond2);

y1Sol = sol.y1;
y2Sol = sol.y2;

disp('Symbolic solution for y1(t):');
disp(y1Sol);
disp('Symbolic solution for y2(t):');
disp(y2Sol);

fplot(y1Sol, [0, 5], 'LineWidth', 2); 
hold on;
fplot(y2Sol, [0, 5], 'LineWidth', 2); 

xlabel('t');
ylabel('Solution Value');
title('Solution of the System Using dsolve');
legend('y_1(t)', 'y_2(t)', 'Location', 'best');
grid on;


%% Exercise 5
% Task 2
clearvars
clc
close all

% a
% Part (a): 

ode_a = @(t, y) [
    (-16/200 * y(1) + 4/200 * y(2));  % dy1/dt
    (16/200 * y(1) - 16/200 * y(2)); % dy2/dt
    % (4 * (y(1)/200 - y(2)/200) + 8 * (0 - y(2)/200))   % dy2/dt
];

y0_a = [10; 20]; 

[t_a, y_a] = ode45(ode_a, [0, 1000], y0_a);

figure;
plot(t_a, y_a(:,1), 'b', 'LineWidth', 2);
hold on;
plot(t_a, y_a(:,2), 'r', 'LineWidth', 2);
xlabel('Time (min)');
ylabel('Salt Amount (kg)');
title('Salt Amount in Tanks (Part a)');
legend('Tank 1', 'Tank 2');
grid on;

%% Part (b):

ode_b = @(t, y) [
    0.5 - 16*y(1)/200 + 4*y(2)/200; % dy1/dt
    16*y(1)/200 - 16*y(2)/200; % dy2/dt
    % (12 * (0 - y(1)/200) + 4 * (y(2)/200 - y(1)/200) + 0.5); 
    % (4 * (y(1)/200 - y(2)/200) + 8 * (0 - y(2)/200))          
];

y0_b = [10; 20];  

[t_b, y_b] = ode45(ode_b, [0, 1000], y0_b);

figure;
plot(t_b, y_b(:,1), 'b', 'LineWidth', 2);
hold on;
plot(t_b, y_b(:,2), 'r', 'LineWidth', 2);
xlabel('Time (min)');
ylabel('Salt Amount (kg)');
title('Salt Amount in Tanks (Part b)');
legend('Tank 1', 'Tank 2');
grid on;

%% Part (d): 

ode_d = @(t, y) [
    (0.02 + 0.01*t - 21*y(1)/200 + 9*y(2)/200);    % dy1/dt
    (16*y(1)/200 - 21*y(2)/200 + 5*y(3)/100);   % dy2/dt
    (5*y(1)/200 - 5*y(3)/100);                 % dy3/dt
    % (12 * (0 - y(1)/200) + 4 * (y(2)/200 - y(1)/200) + 5 * (y(3)/100 - y(1)/200) + 0.5); % dy1/dt
    % (4 * (y(1)/200 - y(2)/200) + 8 * (0 - y(2)/200) + 5 * (y(3)/100 - y(2)/200));       % dy2/dt
    % (5 * (y(1)/200 - y(3)/100) + 5 * (y(2)/200 - y(3)/100))                             % dy3/dt
];

y0_d = [10; 20; 0]; 

[t_d, y_d] = ode45(ode_d, [0, 100], y0_d);

figure;
plot(t_d, y_d(:,1), 'b', 'LineWidth', 2);
hold on;
plot(t_d, y_d(:,2), 'r', 'LineWidth', 2);
plot(t_d, y_d(:,3), 'g', 'LineWidth', 2);
xlabel('Time (min)');
ylabel('Salt Amount (kg)');
title('Salt Amount in Tanks (Part d)');
legend('Tank 1', 'Tank 2', 'Tank 3');
grid on;


%% b
clc; clear; close all;

% (i)
A = [1 0;
     3 4];

[V, D] = eig(A);

% V is the matrix whose columns are eigenvectors
% D is the diagonal matrix of eigenvalues

% (ii) 
d1 = D(1,1);       % first eigenvalue
v1 = V(:,1);       % first eigenvector (column 1 of V)

% (iii) 
v1_length = norm(v1);

% Display the results
disp('Eigenvalues (diagonal of D):');
disp(diag(D));

disp('Eigenvectors (columns of V):');
disp(V);

disp('Extracted eigenvalue d1:');
disp(d1);

disp('Extracted eigenvector v1:');
disp(v1);

disp('Norm of the extracted eigenvector v1:');
disp(v1_length);

% (iv) 
% MATLAB’s eig function normalizes the eigenvectors in V by default,
% but the exact sign or scaling is not unique. Any non-zero scalar 
% multiple of an eigenvector is also a valid eigenvector. MATLAB 
% typically chooses the convention that makes the first element 
% (or a particular element) real and positive, but this is not 
% guaranteed in all cases. What matters is that each column of V 
% is a valid eigenvector corresponding to the eigenvalue in D.


%% Exercise 5
% Task 3
clearvars
clc
close all


% a
m = 1;
k = 0.5;
g = 9.81;

tspan = [0, 10];  

y0 = [100; 0];

f_noDrag = @(t,y) [ y(2); -g ];           

f_withDrag = @(t,y) [ y(2); -g - (k/m)*y(2) ];

[t_noDrag, y_noDrag]     = ode45(f_noDrag,    tspan, y0);
[t_withDrag, y_withDrag] = ode45(f_withDrag,  tspan, y0);

figure('Name','Falling Object - Velocity vs Time','NumberTitle','off');
hold on; grid on;
plot(t_noDrag, y_noDrag(:,2), 'b-', 'LineWidth',2, 'DisplayName','No Drag');
plot(t_withDrag, y_withDrag(:,2), 'r-', 'LineWidth',2, 'DisplayName','With Drag');
xlabel('Time (s)');
ylabel('Velocity v(t) (m/s)');
title('Falling Object: Velocity vs Time');
legend('Location','Best');
hold off;

%% b
m = 1;
k = 0.5;
g = 9.81;

tspan = [0, 10];  
y0 = [100; 0];

f_noDrag = @(t,y) [ y(2); -g ];        
f_withDrag = @(t,y) [ y(2); -g - (k/m)*y(2) ];

[t_noDrag, y_noDrag]     = ode45(f_noDrag,    tspan, y0);
[t_withDrag, y_withDrag] = ode45(f_withDrag,  tspan, y0);

figure('Name','Falling Objects with Velocity Vectors','NumberTitle','off');
hold on; grid on;
axis([tspan(1) tspan(2) 0 110]);
xlabel('Time (s)');
ylabel('Height (m)');
title('Animation: Two Falling Objects with Velocity Arrows (Downward)');

plot(t_noDrag, y_noDrag(:,1), 'b--','HandleVisibility','off');
plot(t_withDrag, y_withDrag(:,1), 'r--','HandleVisibility','off');

markerNoDrag   = plot(t_noDrag(1), y_noDrag(1,1), 'bo','MarkerFaceColor','b','MarkerSize',8, 'DisplayName','No Drag');
markerWithDrag = plot(t_withDrag(1), y_withDrag(1,1), 'ro','MarkerFaceColor','r','MarkerSize',8, 'DisplayName','With Drag');

scaleV = 0.5;  

v0_noDrag   = y_noDrag(1,2);
v0_withDrag = y_withDrag(1,2);
velArrowNoDrag   = quiver(t_noDrag(1), y_noDrag(1,1), 0, v0_noDrag*scaleV, ...
                           'b', 'LineWidth',2, 'MaxHeadSize',0.5, 'AutoScale','off');
velArrowWithDrag = quiver(t_withDrag(1), y_withDrag(1,1), 0, v0_withDrag*scaleV, ...
                           'r', 'LineWidth',2, 'MaxHeadSize',0.5, 'AutoScale','off');


legend('No Drag','With Drag','Location','Best');

tCurrent = 0;     
tEnd     = max(t_noDrag(end), t_withDrag(end));
tStep    = 0.01;  

while tCurrent <= tEnd
    sNoDrag   = interp1(t_noDrag,   y_noDrag(:,1),   tCurrent, 'linear','extrap');
    sWithDrag = interp1(t_withDrag, y_withDrag(:,1), tCurrent, 'linear','extrap');
    
    vNoDrag   = interp1(t_noDrag,   y_noDrag(:,2),   tCurrent, 'linear','extrap');
    vWithDrag = interp1(t_withDrag, y_withDrag(:,2), tCurrent, 'linear','extrap');
    
    set(markerNoDrag,   'XData', tCurrent, 'YData', sNoDrag);
    set(markerWithDrag, 'XData', tCurrent, 'YData', sWithDrag);
    
    set(velArrowNoDrag,   'XData', tCurrent, 'YData', sNoDrag,   'UData', 0, 'VData', vNoDrag*scaleV);
    set(velArrowWithDrag, 'XData', tCurrent, 'YData', sWithDrag, 'UData', 0, 'VData', vWithDrag*scaleV);
    
    pause(tStep);
    tCurrent = tCurrent + tStep;
end

hold off;

%% c
clc; close all; clearvars;

m = 1;
k = 0.5;
g = 9.81;
tspan = [0, 40];  
y0 = [200; 0];

f_noDrag   = @(t,y) [ y(2); -g ];           
f_withDrag = @(t,y) [ y(2); -g - (k/m)*y(2) ];

[t_noDrag, y_noDrag]     = ode45(f_noDrag,   tspan, y0);
[t_withDrag, y_withDrag] = ode45(f_withDrag, tspan, y0);

figure('Name','Falling Objects with Force Vectors','NumberTitle','off');
hold on; grid on;
axis([tspan(1) tspan(2) y0(2) y0(1)]);
xlabel('Time (s)');
ylabel('Height (m)');
title('Falling Objects with Force Vectors');

plot(t_noDrag, y_noDrag(:,1), 'b--','HandleVisibility','off');
plot(t_withDrag, y_withDrag(:,1), 'r--','HandleVisibility','off');

markerNoDrag   = plot(t_noDrag(1),   y_noDrag(1,1), 'bo','MarkerFaceColor','b','MarkerSize',8);
markerWithDrag = plot(t_withDrag(1), y_withDrag(1,1), 'ro','MarkerFaceColor','r','MarkerSize',8);

forceScale = 1;  

F_grav = -m*g;  

offset_noDrag   = 0.0;  
offset_grav_with = 0.0;   
offset_drag_with = 0.0;   
offset_net_with  = 0.0;   


arrowFg_noDrag = quiver(t_noDrag(1)+offset_noDrag, y_noDrag(1,1), 0, F_grav*forceScale, ...
    'b', 'LineWidth',2, 'MaxHeadSize',1, 'AutoScale','off');

F_grav_with = F_grav;  
v0 = y_withDrag(1,2);
F_drag0 = -k * v0;
F_net0  = F_grav_with + F_drag0;

arrowGrav_with = quiver(t_withDrag(1)+offset_grav_with, y_withDrag(1,1), 0, F_grav_with*forceScale, ...
    'r', 'LineWidth',2, 'MaxHeadSize',1, 'AutoScale','off');
arrowDrag_with = quiver(t_withDrag(1)+offset_drag_with, y_withDrag(1,1), 0, F_drag0*forceScale, ...
    'g', 'LineWidth',2, 'MaxHeadSize',1, 'AutoScale','off');
arrowNet_with  = quiver(t_withDrag(1)+offset_net_with, y_withDrag(1,1), 0, F_net0*forceScale, ...
    'k', 'LineWidth',2, 'MaxHeadSize',1, 'AutoScale','off');

tCurrent = 0;     
tEnd     = max(t_noDrag(end), t_withDrag(end));
tStep    = 0.1;  

while tCurrent <= tEnd
    sNoDrag   = interp1(t_noDrag,   y_noDrag(:,1),   tCurrent, 'linear','extrap');
    sWithDrag = interp1(t_withDrag, y_withDrag(:,1), tCurrent, 'linear','extrap');
    
    set(markerNoDrag,   'XData', tCurrent, 'YData', sNoDrag);
    set(markerWithDrag, 'XData', tCurrent, 'YData', sWithDrag);
    
    set(arrowFg_noDrag, 'XData', tCurrent + offset_noDrag, 'YData', sNoDrag, ...
        'UData', 0, 'VData', F_grav*forceScale);
    
    v_with = interp1(t_withDrag, y_withDrag(:,2), tCurrent, 'linear','extrap');
    F_drag = -k * v_with;         
    F_net  = F_grav_with + F_drag;  
    
    set(arrowGrav_with, 'XData', tCurrent + offset_grav_with, 'YData', sWithDrag, ...
        'UData', 0, 'VData', F_grav_with*forceScale);
    set(arrowDrag_with, 'XData', tCurrent + offset_drag_with, 'YData', sWithDrag, ...
        'UData', 0, 'VData', F_drag*forceScale);
    set(arrowNet_with,  'XData', tCurrent + offset_net_with,  'YData', sWithDrag, ...
        'UData', 0, 'VData', F_net*forceScale);
    
    pause(tStep);
    tCurrent = tCurrent + tStep;
end

hold off;
