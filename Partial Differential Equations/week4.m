%% MATLAB 1
clc
close all
clear

% Define two functions C^k(R^2)
syms x y
f = x^2*y + sin(x);
g = exp(x*y);

% Define the order of differentiation for each variable
alpha = [1 2];   % d^3 / (dx dy^2)

% Evaluation point
x_diff = [1 1];

% Compute the derivative using General Leibniz Rule
D_alpha_fg = general_leibniz_rule(f, g, alpha, x_diff, [x,y]);

disp('Result of General Leibniz Rule:')
disp(D_alpha_fg)


% Compute the derivative using General Leibniz Rule
D_alpha_fg = general_leibniz_rule(f, g, alpha, x_diff, [x,y]);

function D_alpha_fg = general_leibniz_rule(f, g, alpha, x_diff, x)
    % Initialize the result
    D_alpha_fg = 0;

    % Generate all multi-indices beta such that beta <= alpha
    betas = generate_multi_indices(alpha);

    % Loop over all beta multi-indices
    for i = 1:size(betas, 1)
        beta = betas(i,:);
        alpha_minus_beta = alpha - beta;

        % Compute the multinomial coefficient
        coeff = multinomial_coeff(alpha, beta);

        % Implement symbolic derivation logic:
        Df = f;
        Dg = g;

        for j = 1:length(x)
            Df = diff(Df, x(j), alpha_minus_beta(j));
            Dg = diff(Dg, x(j), beta(j));
        end

        % Evaluate at differentiation point
        Df_val = subs(Df, x, x_diff);
        Dg_val = subs(Dg, x, x_diff);

        % Update the solution value
        D_alpha_fg = D_alpha_fg + coeff * Df_val * Dg_val;
    end

    % Convert symbolic result to numeric
    D_alpha_fg = double(D_alpha_fg);
end

function coeff = multinomial_coeff(alpha, beta)
    coeff = 1;
    for j = 1:length(alpha)
        coeff = coeff * factorial(alpha(j)) / ...
            (factorial(beta(j)) * factorial(alpha(j) - beta(j)));
    end
end

function indices = generate_multi_indices(alpha)
    % Generate all multi-indices beta such that beta <= alpha
    n = length(alpha);
    indices = [];
    for i = 0:prod(alpha+1)-1      % prod returns product of the elements in its input
        % dec2base(i,j,n) returns number i in base j in minimum n number represenation form
        beta = dec2base(i, max(alpha+1), n) - '0';       % - '0' changes the output from string to number format
        if beta <= alpha
            indices = [indices; beta];
        end
    end
end

%% MATLAB 2
clc
close all
clear

% Define the function
syms x y
f = 0.2*exp(0.5*x + 3*y);

expansion_point = [1 1];
expansion_order = 5;

% Create grid for plotting
[X,Y] = meshgrid(0:0.1:2, 0:0.1:2);

% Plot the original function
F = matlabFunction(f);
figure
surf(X,Y,F(X,Y))
hold on
shading interp
title('Taylor polynomials of different orders')
xlabel('x')
ylabel('y')
zlabel('z')

% Go through each of the different expansion orders and plot each surface
for k = 1:expansion_order
    ts = taylor_series(f, [x,y], expansion_point, k);
    Ts = matlabFunction(ts);
    surf(X,Y,Ts(X,Y))
end

legend('Original f', 'k=1', 'k=2', 'k=3', 'k=4', 'k=5')
hold off



function ts = taylor_series(func, vars, point, k)
    % Initialize output
    taylor_expansion = 0;

    % Generate all values of alphas
    alphas = generate_combs(k);

    % Generate the series
    for i = 1:size(alphas,1)   % Go through each tuple alpha
        cur_alpha = alphas(i,:);

        % Implement derivation logic:
        df = func;
        for j = 1:length(vars)
            df = diff(df, vars(j), cur_alpha(j));
        end

        % Evaluate the value of the derivative at expansion point
        dval = subs(df, vars, point);

        % Generate one term of the Taylor expansion
        alpha_fact = prod(factorial(cur_alpha));
        term = dval / alpha_fact;
        for j = 1:length(vars)
            term = term * (vars(j) - point(j))^cur_alpha(j);
        end

        taylor_expansion = taylor_expansion + term;
    end
    ts = simplify(taylor_expansion);
end



function derivatives = generate_combs(k)
    derivatives = [];
    for alpha1 = 0:k
        for alpha2 = 0:k-alpha1
            derivatives = [derivatives; alpha1, alpha2];
        end
    end
end