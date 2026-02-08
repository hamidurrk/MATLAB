%% MATLAB 1

f = [0.32; 0.25; 0.05];
Aeq = [1 1 1; 0.055 0.035 0];
beq = [100; 4];
lb = [0; 0; 0];
[x, fval] = linprog(f, [], [], Aeq, beq, lb);
x
fval


%% MATLAB 2

% a
f = @(x,y) 1./((2*x-sqrt(2)).^2 + (5*y-pi).^2 + 1);

N = 200;
[j,k] = ndgrid(1:N,1:N);
x = j/N; y = k/N;

M = f(x,y);
heatmap(M);

%b
[~,idx] = max(M(:));
[jN,kN] = ind2sub([N N],idx);

xN = jN/N; 
yN = kN/N;

%c
f = @(x,y) 1./((2*x-sqrt(2)).^2 + (5*y-pi).^2 + 1);

Ns = [20 50 100 200 400];
pts = zeros(numel(Ns),2);

for t = 1:numel(Ns)
    N = Ns(t);
    [j,k] = ndgrid(1:N,1:N);
    M = f(j/N,k/N);
    [~,idx] = max(M(:));
    [jN,kN] = ind2sub([N N],idx);
    pts(t,:) = [jN/N, kN/N];
end

pts

% d
% If the function depends on m variables and each variable is discretized 
% into N points, then the total number of function evaluations is N^M.
% For: m=30; N=1000; we must evaluate 10^90 points, which is more than the
% amount of particles in the observable universe (3.28*10^80).
% Even if every particle in the universe evaluated one point, it would 
% still be insufficient. This nature of dimensionality makes brute-force 
% grid search infeasible for practical high-dimensional optimization 
% problems.