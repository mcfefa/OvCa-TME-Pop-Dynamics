%%%Figuring out how fibroblasts interact with sensitive and resistant cells
%%%MODEL A: F linear term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/TykwithTCAF/results';

coculture_files = {
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/data/Fibro_75S25RTyk_cellsinwell.csv'
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/data/Fibro_50S50RTyk_cellsinwell.csv'
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/data/Fibro_25S75RTyk_cellsinwell.csv'
};

coculture_data = cell(3,1);
S_raw_all = cell(3,1);
R_raw_all = cell(3,1);
ratio_labels = {'75/25 S/R', '50/50 S/R', '25/75 S/R'};

for i = 1:3
    coculture_data{i} = readtable(coculture_files{i});
    S_raw_all{i} = [coculture_data{i}.S_Rep1, coculture_data{i}.S_Rep2, coculture_data{i}.S_Rep3];
    R_raw_all{i} = [coculture_data{i}.R_Rep1, coculture_data{i}.R_Rep2, coculture_data{i}.R_Rep3];
end

t_coculture = coculture_data{1}.Day;
idx_0_10 = (t_coculture >= 0) & (t_coculture <= 10);
t_coculture_10 = t_coculture(idx_0_10);

S_raw_10_all = cell(3,1);
R_raw_10_all = cell(3,1);
for i = 1:3
    S_raw_10_all{i} = S_raw_all{i}(idx_0_10, :);
    R_raw_10_all{i} = R_raw_all{i}(idx_0_10, :);
end

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

%fit for betaS and betaR only before the subtractive or additive term
p0 = [0.5, 0.1];
lb = [0.0, 0.0]; %allow to alternate between positive and negative 
ub = [2.0, 2.0];

res_fun_3pop = @(p) coculture_residuals_3pop_untreated_model1(p, rS, K, alphaS, IC50S, rR, alphaR, IC50R, rF, KF, t_coculture_10, S_raw_10_all, R_raw_10_all);

problem_3pop = createOptimProblem('lsqnonlin', ...
    'x0', p0, ...
    'objective', res_fun_3pop, ...
    'lb', lb, ...
    'ub', ub);

ms_3pop = MultiStart('UseParallel', true, 'Display', 'iter');
nStarts = 200;

[p_hat, resnorm_hat] = run(ms_3pop, problem_3pop, nStarts);

fprintf('Fitted parameters:\n');
fprintf('betaS = %.4f\n', p_hat(1));
fprintf('betaR = %.4f\n', p_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_hat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted parameters: for + for sensitive and + for resistant 
% Fitted parameters:
% betaS = 0.2438
% betaR = 0.5314
% Sum of squared residuals = 810619484102.52

% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted parameters: for - for sensitive and + for resistant 
% betaS = 0.0000
% betaR = 0.2580
% Sum of squared residuals = 1019915864418.05
% IdleTimeout has been reached.
% Parallel pool using the 'Processes' profile is shutting down.

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_log_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_3pop, p_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_log_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n');
fprintf('betaS: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('betaR: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
% 95 Confidence Intervals:
% betaS: [0.1717, 0.3159]
% betaR: [0.4343, 0.6285]


[~, BIC_conditions] = coculture_residuals_3pop_untreated_model1(p_hat, rS, ...
    K, alphaS, IC50S, rR, alphaR, IC50R, rF, KF, ...
    t_coculture_10, S_raw_10_all, R_raw_10_all);
ratio_labels = {'75/25 S/R', '50/50 S/R', '25/75 S/R'};

fprintf('\nBIC per condition:\n');
for i = 1:length(BIC_conditions)
    fprintf('  %s: BIC = %.2f\n', ratio_labels{i}, BIC_conditions(i));
end
% BIC per condition:
%   75/25 S/R: BIC = 1490.94
%   50/50 S/R: BIC = 1463.04
%   25/75 S/R: BIC = 1445.31

figure;
for i = 1:3
    S_raw = S_raw_10_all{i};
    R_raw = R_raw_10_all{i};
    [nTimes, nReps] = size(S_raw);

    S_mod_all = zeros(nTimes, nReps);
    R_mod_all = zeros(nTimes, nReps);

    for j = 1:nReps
        S0 = S_raw(1,j);
        R0 = R_raw(1,j);
        F0 = 10000;

        [S_mod, R_mod, F_mod] = simulate_coculture_3pop_untreated_model1(t_coculture_10, S0, R0, F0, p_hat, rS, ...
    K, alphaS, IC50S, rR, alphaR, IC50R, rF, KF);
        S_mod_all(:,j) = S_mod;
        R_mod_all(:,j) = R_mod;
    end

    S_mod_mean = mean(S_mod_all, 2);
    R_mod_mean = mean(R_mod_all, 2);

    subplot(1,3,i); hold on;
    plot(t_coculture_10, S_raw, 'ro', 'MarkerFaceColor','r', 'MarkerSize',4);
    plot(t_coculture_10, R_raw, 'go', 'MarkerFaceColor','g', 'MarkerSize',4);
    plot(t_coculture_10, S_mod_mean, 'r-', 'LineWidth',2);
    plot(t_coculture_10, R_mod_mean, 'g-', 'LineWidth',2);

    title(ratio_labels{i});
    xlabel('Time (days)');
    ylabel('Cell count');
    ylim([0 1400000]);
    grid on;
    legend({'S data','R data','S model','R model'}, 'Location','best');
end

sgtitle(sprintf('3-population untreated fit'));

saveas(gcf, fullfile(outdir,'co-culture_3pop_untreated_fit_model1.png'));


figure;
for i = 1:3
    S_raw = S_raw_10_all{i};
    R_raw = R_raw_10_all{i};
    [nTimes, nReps] = size(S_raw);

    S_mod_all = zeros(nTimes, nReps);
    R_mod_all = zeros(nTimes, nReps);
    F_mod_all = zeros(nTimes, nReps);

    for j = 1:nReps
        S0 = S_raw(1,j);
        R0 = R_raw(1,j);
        F0 = 10000;

        [S_mod, R_mod, F_mod] = simulate_coculture_3pop_untreated_model1( ...
            t_coculture_10, S0, R0, F0, p_hat, rS, K, alphaS, IC50S, ...
            rR, alphaR, IC50R, rF, KF);

        S_mod_all(:,j) = S_mod;
        R_mod_all(:,j) = R_mod;
        F_mod_all(:,j) = F_mod;
    end

    S_mean = mean(S_raw, 2);
    S_sd   = std(S_raw, 0, 2);
    R_mean = mean(R_raw, 2);
    R_sd   = std(R_raw, 0, 2);

    S_mod_mean = mean(S_mod_all, 2);
    R_mod_mean = mean(R_mod_all, 2);
    F_mod_mean = mean(F_mod_all, 2);

    Tot_mean = S_mean + R_mean;
    Tot_mod_mean = S_mod_mean + R_mod_mean + F_mod_mean;

    writematrix([t_coculture_10, S_mod_all], ...
        fullfile(outdir, sprintf('S_sim_i%d_untreated_3pop_10days_model1.csv', i)));

    writematrix([t_coculture_10, R_mod_all], ...
        fullfile(outdir, sprintf('R_sim_i%d_untreated_3pop_10days_model1.csv', i)));

    writematrix([t_coculture_10, F_mod_all], ...
        fullfile(outdir, sprintf('F_sim_i%d_untreated_3pop_10days_model1.csv', i)));

    subplot(2,2,i); hold on;

    Sdata = errorbar(t_coculture_10, S_mean, S_sd, 'r.', 'MarkerSize', 12, 'LineStyle','none');
    Rdata = errorbar(t_coculture_10, R_mean, R_sd, 'g.', 'MarkerSize', 12, 'LineStyle','none');
    Tdata = plot(t_coculture_10, Tot_mean, 'b.', 'MarkerSize', 12);

    Smodel = plot(t_coculture_10, S_mod_mean, 'r-', 'LineWidth', 2);
    Rmodel = plot(t_coculture_10, R_mod_mean, 'g-', 'LineWidth', 2);
    Fmodel = plot(t_coculture_10, F_mod_mean, 'm-', 'LineWidth', 2);
    Tmodel = plot(t_coculture_10, Tot_mod_mean, 'b-', 'LineWidth', 2);

    title(sprintf('%s', ratio_labels{i}));
    xlabel('Time (days)');
    ylabel('Cell count');
    ylim([0 1400000]);
    grid on;
    legend([Sdata(1), Rdata(1), Tdata(1), Smodel(1), Rmodel(1), Fmodel(1), Tmodel(1)], ...
        {'S data','R data','Total data','S model','R model','F model','Total model'}, ...
        'Location','best');
end

sgtitle('Untreated 3-population fit');

saveas(gcf, fullfile(outdir,'co-culture_3pop_untreated_fitwithallpop_model1.png'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%MODEL B: F affects carrying capacity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/TykwithTCAF/results';

coculture_files = {
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/data/Fibro_75S25RTyk_cellsinwell.csv'
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/data/Fibro_50S50RTyk_cellsinwell.csv'
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/data/Fibro_25S75RTyk_cellsinwell.csv'
};

coculture_data = cell(3,1);
S_raw_all = cell(3,1);
R_raw_all = cell(3,1);
ratio_labels = {'75/25 S/R', '50/50 S/R', '25/75 S/R'};

for i = 1:3
    coculture_data{i} = readtable(coculture_files{i});
    S_raw_all{i} = [coculture_data{i}.S_Rep1, coculture_data{i}.S_Rep2, coculture_data{i}.S_Rep3];
    R_raw_all{i} = [coculture_data{i}.R_Rep1, coculture_data{i}.R_Rep2, coculture_data{i}.R_Rep3];
end

t_coculture = coculture_data{1}.Day;
idx_0_10 = (t_coculture >= 0) & (t_coculture <= 10);
t_coculture_10 = t_coculture(idx_0_10);

S_raw_10_all = cell(3,1);
R_raw_10_all = cell(3,1);
for i = 1:3
    S_raw_10_all{i} = S_raw_all{i}(idx_0_10, :);
    R_raw_10_all{i} = R_raw_all{i}(idx_0_10, :);
end

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

%fit for betaS and betaR only before the subtractive or additive term
p0 = [0.5, 0.1];
lb = [0.0, 0.0]; %allow to alternate between positive and negative 
ub = [5.0, 5.0];

res_fun_3pop = @(p) coculture_residuals_3pop_untreated_model3(p, rS, K, alphaS, IC50S, rR, alphaR, IC50R, rF, KF, t_coculture_10, S_raw_10_all, R_raw_10_all);

problem_3pop = createOptimProblem('lsqnonlin', ...
    'x0', p0, ...
    'objective', res_fun_3pop, ...
    'lb', lb, ...
    'ub', ub);

ms_3pop = MultiStart('UseParallel', true, 'Display', 'iter');
nStarts = 200;

[p_hat, resnorm_hat] = run(ms_3pop, problem_3pop, nStarts);

fprintf('Fitted parameters:\n');
fprintf('betaS = %.4f\n', p_hat(1));
fprintf('betaR = %.4f\n', p_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_hat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted parameters: + sensitive and - for resistant 
% betaS = 0.0000
% betaR = 1.3182
% Sum of squared residuals = 711897456791.07
% IdleTimeout has been reached.

% 197 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted parameters: - sensitive and - resistant
% betaS = 1.3613
% betaR = 2.0631
% Sum of squared residuals = 371529973794.25

% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted parameters: for + sensitive and + resistant 
% betaS = 0.4375
% betaR = 0.0000
% Sum of squared residuals = 1648313582107.38


% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_log_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_3pop, p_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_log_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n');
fprintf('betaS: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('betaR: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
% 95 Confidence Intervals:

[~, BIC_conditions] = coculture_residuals_3pop_untreated_model3(p_hat, rS, ...
    K, alphaS, IC50S, rR, alphaR, IC50R, rF, KF, ...
    t_coculture_10, S_raw_10_all, R_raw_10_all);
ratio_labels = {'75/25 S/R', '50/50 S/R', '25/75 S/R'};

fprintf('\nBIC per condition:\n');
for i = 1:length(BIC_conditions)
    fprintf('  %s: BIC = %.2f\n', ratio_labels{i}, BIC_conditions(i));
end
% BIC per condition: + for sensitive and - for resistant 
%   75/25 S/R: BIC = 1468.06
%   50/50 S/R: BIC = 1479.03
%   25/75 S/R: BIC = 1421.98

% BIC per condition: - sensitive and - resistant 
%   75/25 S/R: BIC = 1426.89
%   50/50 S/R: BIC = 1417.99
%   25/75 S/R: BIC = 1406.51

figure;
for i = 1:3
    S_raw = S_raw_10_all{i};
    R_raw = R_raw_10_all{i};
    [nTimes, nReps] = size(S_raw);

    S_mod_all = zeros(nTimes, nReps);
    R_mod_all = zeros(nTimes, nReps);

    for j = 1:nReps
        S0 = S_raw(1,j);
        R0 = R_raw(1,j);
        F0 = 10000;

        [S_mod, R_mod, F_mod] = simulate_coculture_3pop_untreated_model3(t_coculture_10, S0, R0, F0, p_hat, rS, ...
    K, alphaS, IC50S, rR, alphaR, IC50R, rF, KF);
        S_mod_all(:,j) = S_mod;
        R_mod_all(:,j) = R_mod;
    end

    S_mod_mean = mean(S_mod_all, 2);
    R_mod_mean = mean(R_mod_all, 2);

    subplot(1,3,i); hold on;
    plot(t_coculture_10, S_raw, 'ro', 'MarkerFaceColor','r', 'MarkerSize',4);
    plot(t_coculture_10, R_raw, 'go', 'MarkerFaceColor','g', 'MarkerSize',4);
    plot(t_coculture_10, S_mod_mean, 'r-', 'LineWidth',2);
    plot(t_coculture_10, R_mod_mean, 'g-', 'LineWidth',2);

    title(ratio_labels{i});
    xlabel('Time (days)');
    ylabel('Cell count');
    ylim([0 1400000]);
    grid on;
    legend({'S data','R data','S model','R model'}, 'Location','best');
end

sgtitle(sprintf('3-population untreated fit'));

saveas(gcf, fullfile(outdir,'co-culture_3pop_untreated_fit_model3.png'));


figure;
for i = 1:3
    S_raw = S_raw_10_all{i};
    R_raw = R_raw_10_all{i};
    [nTimes, nReps] = size(S_raw);

    S_mod_all = zeros(nTimes, nReps);
    R_mod_all = zeros(nTimes, nReps);
    F_mod_all = zeros(nTimes, nReps);

    for j = 1:nReps
        S0 = S_raw(1,j);
        R0 = R_raw(1,j);
        F0 = 10000;

        [S_mod, R_mod, F_mod] = simulate_coculture_3pop_untreated_model3( ...
            t_coculture_10, S0, R0, F0, p_hat, rS, K, alphaS, IC50S, ...
            rR, alphaR, IC50R, rF, KF);

        S_mod_all(:,j) = S_mod;
        R_mod_all(:,j) = R_mod;
        F_mod_all(:,j) = F_mod;
    end

    S_mean = mean(S_raw, 2);
    S_sd   = std(S_raw, 0, 2);
    R_mean = mean(R_raw, 2);
    R_sd   = std(R_raw, 0, 2);

    S_mod_mean = mean(S_mod_all, 2);
    R_mod_mean = mean(R_mod_all, 2);
    F_mod_mean = mean(F_mod_all, 2);

    Tot_mean = S_mean + R_mean;
    Tot_mod_mean = S_mod_mean + R_mod_mean + F_mod_mean;

    writematrix([t_coculture_10, S_mod_all], ...
        fullfile(outdir, sprintf('S_sim_i%d_untreated_3pop_10days_model3.csv', i)));

    writematrix([t_coculture_10, R_mod_all], ...
        fullfile(outdir, sprintf('R_sim_i%d_untreated_3pop_10days_model3.csv', i)));

    writematrix([t_coculture_10, F_mod_all], ...
        fullfile(outdir, sprintf('F_sim_i%d_untreated_3pop_10days_model3.csv', i)));

    subplot(2,2,i); hold on;

    Sdata = errorbar(t_coculture_10, S_mean, S_sd, 'r.', 'MarkerSize', 12, 'LineStyle','none');
    Rdata = errorbar(t_coculture_10, R_mean, R_sd, 'g.', 'MarkerSize', 12, 'LineStyle','none');
    Tdata = plot(t_coculture_10, Tot_mean, 'b.', 'MarkerSize', 12);

    Smodel = plot(t_coculture_10, S_mod_mean, 'r-', 'LineWidth', 2);
    Rmodel = plot(t_coculture_10, R_mod_mean, 'g-', 'LineWidth', 2);
    Fmodel = plot(t_coculture_10, F_mod_mean, 'm-', 'LineWidth', 2);
    Tmodel = plot(t_coculture_10, Tot_mod_mean, 'b-', 'LineWidth', 2);

    title(sprintf('%s', ratio_labels{i}));
    xlabel('Time (days)');
    ylabel('Cell count');
    ylim([0 1600000]);
    grid on;
    legend([Sdata(1), Rdata(1), Tdata(1), Smodel(1), Rmodel(1), Fmodel(1), Tmodel(1)], ...
        {'S data','R data','Total data','S model','R model','F model','Total model'}, ...
        'Location','best');
end

sgtitle('Untreated 3-population fit');

saveas(gcf, fullfile(outdir,'co-culture_3pop_untreated_fitwithallpop_model3.png'));
