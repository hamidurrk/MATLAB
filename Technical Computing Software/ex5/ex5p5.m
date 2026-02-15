function final_conc = ex5p5()

    k = 0.02;           % Rate constant (min^-1)
    f0 = 20;            % Initial concentration (mg/mL)
    t_final = 100;      % Total time (min)
    rel_tol = 0.05;     % 5% relative accuracy

    final_conc = euler_refine(@(f) -k*f, f0, t_final, rel_tol);

end

function f_end = euler_refine(dfdt, f0, t_final, rel_tol)
    N = 10; % initial number of steps
    prev = 0;
    while true
        h = t_final / N;
        f = f0;
        for i = 1:N
            f = f + h * dfdt(f);
        end
        if N > 10
            if abs(f - prev)/abs(f) < rel_tol
                break;
            end
        end
        prev = f;
        N = N * 2;
    end
    f_end = f;
end