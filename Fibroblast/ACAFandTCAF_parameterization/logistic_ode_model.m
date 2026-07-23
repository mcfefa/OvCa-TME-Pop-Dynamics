function y_model = logistic_ode_model(p, t_days, y_raw)
    % p = [K, r]
    K = p(1);
    r = p(2);

    [nDays, nReps] = size(y_raw);
    t_all = repmat(t_days, nReps, 1);

    y_model = zeros(numel(t_all), 1);
    idx = 1;

    for j = 1:nReps
        % Initial condition for this replicate = observation at day 0
        C0 = y_raw(1, j);

        % Solve ODE over days 1–14
        % Use continuous time from 1 to 14 and then pick values at each day
        tspan = [t_days(1), t_days(end)];

        ode_fun = @(t, C) r * C * (1 - C / K);
        opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t_sol, C_sol] = ode45(ode_fun, tspan, C0, opts);

        % Interpolate at integer days (0–14)
        C_at_days = interp1(t_sol, C_sol, t_days, 'pchip');

        % Fill into y_model
        y_model(idx:idx+nDays-1) = C_at_days;
        idx = idx + nDays;
    end
end