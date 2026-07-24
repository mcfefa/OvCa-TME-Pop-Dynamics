function [S_model, R_model] = simulate_coculture_3pop_untreated_model2(t_data, S0, R0, p, K, alphaS, nS, IC50S, alphaR, nR, IC50R)

    rS = p(1);
    rR = p(2);

    tspan = [min(t_data), max(t_data)];
    ode_fun = @(t,y) coculture_ode_3pop_untreated_model2(t, y, rS, K, alphaS, nS, IC50S, rR, alphaR, nR, IC50R);
    opts = odeset('RelTol',1e-8, 'AbsTol',1e-10);

    [t_sol, y_sol] = ode15s(ode_fun, tspan, [S0, R0], opts);

    S_model = interp1(t_sol, y_sol(:,1), t_data, 'pchip');
    R_model = interp1(t_sol, y_sol(:,2), t_data, 'pchip');
end