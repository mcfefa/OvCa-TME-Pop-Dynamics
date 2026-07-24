function result = run_coculture_constant_simulation_new(Ce, S0, R0, F0)
%parameters from fitting all data to most sensitive parameters 4

    %parameters from finalized model 
    rS = 0.6066;
    K = 469571.0714;
    alphaS = 0.1713;
    IC50S = 1.974574;
    
    rR = 0.6454;
    alphaR = 0.0006;
    IC50R = 4.210567;
    
    %parameters from fitting the fibroblast data
    rF = 0.4618;
    KF = 163096.7334;
    
    %interactions with fibroblasts
    betaS = 1.3613;
    betaR = 2.0631;
    
    % Time span
    tspan = 0:0.2:25;

    % Initial conditions
    y0 = [S0; R0; F0];

    % ODE function handle
    odefun = @(t,y) coculture_ode_sens2_m13( ...
        t, y, rS, K, alphaS, IC50S, rR, alphaR, IC50R, Ce, betaS, betaR, rF, KF);

    % Solve
    [t, y] = ode15s(odefun, tspan, y0);

    % Extract populations
    S = y(:,1);
    R = y(:,2);
    F = y(:,3);
    T = S + R;

    % Save results to table
    result = table(t, S, R, F, T, ...
        'VariableNames', {'time', 'S', 'R', 'F', 'T'});

    % Create output folder if needed
    if ~exist('output', 'dir')
        mkdir('output');
    end

    % Save CSV
    csvname = sprintf('output/coculture_new_Ce_%g_S0_%g_R0_%g_F0_%g.csv', Ce, S0, R0, F0);
    writetable(result, csvname);

    % Plot
    fig = figure('Color','w','Position',[100 100 900 600]);
    plot(t, S, 'r', 'LineWidth', 2); hold on;
    plot(t, R, 'g', 'LineWidth', 2);
    plot(t, F, 'm', 'LineWidth', 2);
    plot(t, T, 'b', 'LineWidth', 2);
    xlabel('Time');
    ylabel('Cells');
    legend('Sensitive (S)', 'Resistant (R)', 'Fibroblast (F)', 'Total (T)', 'Location', 'best');
    title(sprintf('Co-culture dynamics: Ce = %g, S0 = %g, R0 = %g, F0 = %g', Ce, S0, R0, F0));
    grid on;

    % Save figure
    pngname = sprintf('output/coculture_new_Ce_%g_S0_%g_R0_%g_F0_%g.png', Ce, S0, R0, F0);
    exportgraphics(fig, pngname, 'Resolution', 300);
end

%to run this function do this result = run_coculture_constant_simulation_new(1,22500,7500, 10000)