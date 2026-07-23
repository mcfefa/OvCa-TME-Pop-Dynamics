function y_model = genlog_ode_model(p, t_days, y_raw)
    % p = [K, r, nu]
    K  = p(1);
    r  = p(2);
    nu = p(3);

    [nDays, nReps] = size(y_raw);
    t_all = repmat(t_days, nReps, 1);

    y_model = zeros(numel(t_all), 1);
    idx = 1;

    for j = 1:nReps
        C0 = y_raw(1, j);               % initial value for replicate j

        tspan = [t_days(1), t_days(end)];
        ode_fun = @(t, C) r * C * (1 - (C / K).^nu);
        opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t_sol, C_sol] = ode45(ode_fun, tspan, C0, opts);

        C_at_days = interp1(t_sol, C_sol, t_days, 'pchip');

        y_model(idx:idx+nDays-1) = C_at_days;
        idx = idx + nDays;
    end
end