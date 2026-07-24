function result = run_coculture_adaptive_simulation3_new(S0, R0, F0, u0, threshFrac)
    % Inputs:
    % S0, R0     initial sensitive and resistant cell counts
    % u0         initial dose
    % threshFrac threshold for percent increase (e.g. 0.5 for 50%)

    %parameters from fitting all data to most sensitive parameters
    p.rS = 0.7719;
    p.K = 1112765.1958;
    p.alphaS = 0.0911;
    p.nS = 3.7779;
    p.IC50S = 0.9960177;

    p.rR = 0.6319;
    p.alphaR = 1.7982;
    p.nR = 2.5166;
    p.IC50R = 8.338535;

    %parameters from fitting the fibroblast data
    p.rF = 0.7817;
    p.KF = 109300.8263;
    
    %interactions with fibroblasts
    p.betaS = 0.5688;
    p.betaR = 0.0076;
    
    % Time grid: save every 0.2 days
    tspan = 0:0.2:25;

    % Step schedule: 0, 0.25, ..., 1.5
    dose_levels = 0:0.25:1.5;
    umax = 1.5;
    umin = 0.0;

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
            current_u = umax;
        else
            current_u = 0;
        end
        % Snap to allowed dose levels
        %[~, ix] = min(abs(dose_levels - current_u));
        %current_u = dose_levels(ix);

        last_day_T = day_end_T;
    end

    % Save results
    result = table(tspan(:), S, R, F, T, dose, ...
        'VariableNames', {'time', 'S', 'R', 'F', 'T', 'dose_u1'});

    if ~exist('output_adaptive3', 'dir')
        mkdir('output_adaptive3');
    end

    csvname = sprintf('output_adaptive3/adaptive_new_S0_%0.0f_R0_%0.0f_F0_%0.0f_u0_%0.2f_thresh_%0.2f.csv', ...
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
    ylim([umin umax]);
    yticks(dose_levels);
    title('Adaptive dose schedule');
    grid on;

    png = sprintf('output_adaptive3/adaptive_new_S0_%0.0f_R0_%0.0f_F0_%0.0f_u0_%0.2f_thresh_%0.2f.png', ...
        S0, R0, F0, u0, threshFrac);
    exportgraphics(fig, png, 'Resolution', 300);
end

%to run the function do this resul =
%run_coculture_adaptive_simulation3_new(22500, 7500, 10000, 0.0, 0.5);

%max 1.0, thres 0.1