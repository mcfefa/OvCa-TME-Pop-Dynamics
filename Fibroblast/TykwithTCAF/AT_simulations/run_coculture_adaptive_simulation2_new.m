function result = run_coculture_adaptive_simulation2_new(S0, R0, F0, u0, threshFrac)
    % Inputs:
    % S0, R0     initial sensitive and resistant cell counts
    % u0         initial dose
    % threshFrac threshold for percent increase (e.g. 0.5 for 50%)

    %parameters from finalized model 
    p.rS = 0.6066;
    p.K = 469571.0714;
    p.alphaS = 0.1713;
    p.IC50S = 1.974574;
    
    p.rR = 0.6454;
    p.alphaR = 0.0006;
    p.IC50R = 4.210567;
    
    %parameters from fitting the fibroblast data
    p.rF = 0.4618;
    p.KF = 163096.7334;
    
    %interactions with fibroblasts
    p.betaS = 1.3613;
    p.betaR = 2.0631;

    % Time grid: save every 0.2 days
    tspan = 0:0.2:25;

    % Step schedule: 0, 0.25, ..., 1.5
    dose_levels = 0:0.25:1.5;

    % Initial conditions
    y0 = [S0; R0; F0];

    % Preallocate
    n = numel(tspan);
    S = zeros(n,1);
    R = zeros(n,1);
    F = zeros(n,1);
    T = zeros(n,1);
    dose = zeros(n,1);

    S(1) = S0;
    R(1) = R0;
    F(1) = F0;
    T(1) = S0 + R0;
    dose(1) = u0;

    current_u = u0;
    last_day_T = T(1);

    % Day-by-day simulation
    for day = 1:25
        idx_start = find(abs(tspan - (day-1)) < 1e-12, 1);
        idx_end   = find(abs(tspan - day) < 1e-12, 1);

        y0_day = [S(idx_start); R(idx_start); F(idx_start)];

        odefun = @(t,y) coculture_ode_adaptive( ...
            t, y, p, current_u);

        [tloc, yloc] = ode15s(odefun, [day-1 day], y0_day);

        % Interpolate onto the fixed 0.2-day grid
        S(idx_start:idx_end) = interp1(tloc, yloc(:,1), tspan(idx_start:idx_end), 'pchip');
        R(idx_start:idx_end) = interp1(tloc, yloc(:,2), tspan(idx_start:idx_end), 'pchip');
        F(idx_start:idx_end) = interp1(tloc, yloc(:,3), tspan(idx_start:idx_end), 'pchip');
        T(idx_start:idx_end) = S(idx_start:idx_end) + R(idx_start:idx_end);
        dose(idx_start:idx_end) = current_u;

        % End-of-day population change
        day_end_T = T(idx_end);
        rel_change = (day_end_T - last_day_T) / last_day_T;

        % Adaptive dose update for next day
        if rel_change >= threshFrac
            current_u = min(1.5, current_u + 0.25);
        elseif rel_change < 0
            current_u = max(0, current_u - 0.25);
        end

        % Snap to allowed dose levels
        [~, ix] = min(abs(dose_levels - current_u));
        current_u = dose_levels(ix);

        last_day_T = day_end_T;
    end

    % Save results
    result = table(tspan(:), S, R, F, T, dose, ...
        'VariableNames', {'time', 'S', 'R', 'F', 'T', 'dose_u1'});

    if ~exist('output_adaptive2', 'dir')
        mkdir('output_adaptive2');
    end

    csvname = sprintf('output_adaptive2/adaptive_new_S0_%0.0f_R0_%0.0f_F0_%0.0f_u0_%0.2f_thresh_%0.2f.csv', ...
        S0, R0, F0, u0, threshFrac);
    writetable(result, csvname);

    % Plot populations
    fig = figure('Color','w','Position',[100 100 900 750]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    nexttile
    plot(tspan, S, 'r', 'LineWidth', 2); hold on;
    plot(tspan, R, 'g', 'LineWidth', 2);
    plot(tspan, F, 'm', 'LineWidth', 2);
    plot(tspan, T, 'b', 'LineWidth', 2);
    xlabel('Time (days)');
    ylabel('Cells');
    legend('Sensitive (S)', 'Resistant (R)', 'Fibroblast (F)', 'Total (T)', 'Location', 'best');
    title(sprintf('Adaptive therapy: S0=%g, R0=%g, F0=%g', S0, R0, F0));
    grid on;



    % Plot dose
    nexttile
    stairs(tspan, dose, 'k', 'LineWidth', 2);
    xlabel('Time (days)');
    ylabel('Dose u_1');
    ylim([0 1.5]);
    yticks(dose_levels);
    title('Adaptive dose schedule');
    grid on;

    png = sprintf('output_adaptive2/adaptive_new_S0_%0.0f_R0_%0.0f_F0_%0.0f_u0_%0.2f_thresh_%0.2f.png', ...
        S0, R0, F0, u0, threshFrac);
    exportgraphics(fig, png, 'Resolution', 300);
end

%to run the function do this resul =
%run_coculture_adaptive_simulation2_new(22500, 7500, 10000, 0.0, 0.5);

