%%Trying to estimate parameters with MATLAB for fibroblast cells (ACAF and
%%TCAF)
%Compare logistic and generalized logistic models 

%Start with ACAF for 10 days 
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/ACAFandTCAF_parameterization/results';
%%%%Sensitive cells first 
ACAF_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/data/ACAF_untreated_40k_cellsinwell.csv');
ACAF_t_days = ACAF_data.Day;
ACAF_10t_days = ACAF_t_days(1:11);
ACAF_y_raw = [ACAF_data.Rep1, ACAF_data.Rep2, ACAF_data.Rep3];
ACAF_10y_raw = ACAF_y_raw(1:11,:);
[nDays, nReps] = size(ACAF_10y_raw);

%Vectorize time and data so all replicates are used, no need to do this to
%use lsqnonlin
ACAF_10t_all = repmat(ACAF_10t_days, nReps, 1); %should be 45 x 1 
ACAF_10y_all = ACAF_10y_raw(:); %should be 45 x 1 
ACAF_10N = numel(ACAF_10y_all); %total observations so 45


%Define objective function for both models, need to define a residual
%function to be able to use lsqnonlin
logistic_res = @(p, t_dummy) logistic_ode_model(p, ACAF_10t_days, ACAF_10y_raw)-ACAF_10y_all;

genlog_res   = @(p, t_dummy) genlog_ode_model(p, ACAF_10t_days, ACAF_10y_raw)-ACAF_10y_all;

%Define initial guesses for logistic
K0 = 1.4e6; %max(y_all); but it could be higher than the last day counts
r0 = 0.5; %pretty small

p0_log = [K0, r0];
lb_log = [0,  0];     % K >= 0, r >= 0
ub_log = [Inf, Inf]; %could define a differen upper bound

% Create lsqcurvefit problem for logistic
problem_log = createOptimProblem( ...
    'lsqnonlin', ... %eve though we are using a multistart, we still need to define the function to minimize
    'x0',        p0_log, ...
    'objective', logistic_res, ...
    'lb',        lb_log, ...
    'ub',        ub_log);

ms = MultiStart('UseParallel', true,'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % set false if no Parallel Toolbox
nStarts = 200; %could do more starts if we wanted to 

[p_log, resnorm_log] = run(ms, problem_log, nStarts); %solve the problem with multistart 
fprintf('Fitted logistic parameters for ACAF (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_log(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_log(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_log);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted logistic parameters for ACAF (10 days)
% Fitted growth rate r = 0.7817
% Fitted carrying capacity k = 109300.8263
% Sum of squared residuals = 31214277739.31

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_log_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(logistic_res, p_log, lb_log, ub_log);

% Compute 95% confidence intervals
ci = nlparci(p_log_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95% Confidence Intervals:\n');
fprintf('r: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
fprintf('K: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% <stopping criteria details>
% 95r: [0.5054, 1.0579]
% K: [89192.3271, 129407.5160]

%Define initial guesses for generalized logistic
K0  = 1.4e6;         %max(y_all);
r0  = 0.5;
nu0 = 1;            % nu=1 reduces to logistic

p0_gen = [K0, r0, nu0];
lb_gen = [0,  0,   0.1];   % nu > 0
ub_gen = [Inf, Inf, 30];

problem_gen = createOptimProblem( ...
    'lsqnonlin', ...
    'x0',        p0_gen, ...
    'objective', genlog_res, ...
    'lb',        lb_gen, ...
    'ub',        ub_gen);

[p_gen, resnorm_gen] = run(ms, problem_gen, nStarts);
fprintf('Fitted generalized logistic parameters for ACAF (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_gen(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_gen(1));
fprintf('Fitted n parameter n = %.4f\n', p_gen(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_gen);
% 199 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted generalized logistic parameters for ACAF (10 days)
% Fitted growth rate r = 0.5938
% Fitted carrying capacity k = 516423.7758
% Fitted n parameter n = 1.4037
% Sum of squared residuals = 17332788874.74

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_gen_final, resnorm_gen_final, residual_gen_final, exitflag_gen, output_gen, lambda_gen, J_gen] = ...
lsqnonlin(genlog_res, p_gen, lb_gen, ub_gen);

% Compute 95% confidence intervals
ci_gen = nlparci(p_gen_final, residual_gen_final, 'jacobian', J_gen);

% Print confidence intervals
fprintf('95% Confidence Intervals:\n');
fprintf('r: [%.4f, %.4f]\n', ci_gen(2,1), ci_gen(2,2));
fprintf('K: [%.4f, %.4f]\n', ci_gen(1,1), ci_gen(1,2));
fprintf('nu: [%.4f, %.4f]\n', ci_gen(3,1), ci_gen(3,2)); 
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% <stopping criteria details>
% 95r: [0.5285, 0.6591] (0.6591-0.5285)/0.5938 = 0.2199
% K: [467927.6554, 564919.8963] (564919.8963-467927.6554)/516423.7758 =
% 0.1878
% nu: [0.6465, 2.1607] (2.1607-0.6465)/1.4037 = 1.0787

%%Plot BIC for the models 
N = numel(ACAF_10y_all);  
k_log = length(p_log);
k_gen = length(p_gen);
SSE_log = resnorm_log; %here the resnorm is already squared
SSE_gen = resnorm_gen;

% BIC formula: N*ln(SSE/N) + k*ln(N)
BIC_log = N * log(SSE_log / N) + k_log * log(N);
BIC_gen = N * log(SSE_gen / N) + k_gen * log(N);

% Display results
fprintf('Logistic:     BIC = %.3f\n', BIC_log); %Logistic:     BIC = 689.025
fprintf('Gen Logistic: BIC = %.3f\n', BIC_gen);

if BIC_log < BIC_gen
    fprintf('Logistic preferred (lower BIC)\n');
else
    fprintf('Gen Logistic preferred (lower BIC)\n');
end

% Logistic:     BIC = 671.077
% Gen Logistic: BIC = 673.108
% Logistic preferred (lower BIC)

%%Plot the curves
% Dense time grid for plotting continuous curves
t_fine = linspace(0, 14, 200)';

% Predict for each replicate separately at t_fine, using fitted params
% Example: use the same ODE solvers but with t_fine instead of 1–14
% (for simplicity, here we interpolate from 1–14)

y_log_fit_full = logistic_ode_model(p_log, ACAF_10t_days, ACAF_10y_raw);
y_gen_fit_full = genlog_ode_model(p_gen, ACAF_10t_days, ACAF_10y_raw);

% Reshape back to 14×3 for plotting
y_log_fit = reshape(y_log_fit_full, [nDays, nReps]);
y_gen_fit = reshape(y_gen_fit_full, [nDays, nReps]);

figure;
hold on;
% Plot data
plot(ACAF_10t_days, ACAF_10y_raw, 'ko', 'MarkerFaceColor',[0.8 0.8 0.8]);

% Plot fitted curves per replicate (logistic)
for j = 1:nReps
    plot(ACAF_10t_days, y_log_fit(:,j), 'b-', 'LineWidth', 1.5);
end

% Plot generalized logistic as dashed red (optional)
for j = 1:nReps
    plot(ACAF_10t_days, y_gen_fit(:,j), 'r--', 'LineWidth', 1.5);
end

xlabel('Time (days)');
ylabel('Cell count (or density)');
legend('Data','Logistic fit','Gen. logistic fit','Location','best');
grid on;
hold off;

saveas(gcf, fullfile(outdir,'ACAF_untreated_logandgenlog_10days.png'));

writematrix(y_log_fit, fullfile(outdir,'ACAF_untreated_log_10days.csv'));
writematrix(y_gen_fit, fullfile(outdir,'ACAF_untreated_genlog_10days.csv'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Repeat the 10 day process for TCAF with lsqnonlin
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/ACAFandTCAF_parameterization/results';

%%%%Resistant cells cells  
TCAF_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Fibroblast/data/TCAF_untreated_40k_cellsinwell.csv');
TCAF_t_days = TCAF_data.Day;
TCAF_10t_days = TCAF_t_days(1:11);
TCAF_y_raw = [TCAF_data.Rep1, TCAF_data.Rep2, TCAF_data.Rep3];
TCAF_10y_raw = TCAF_y_raw(1:11,:);
[nDays, nReps] = size(TCAF_10y_raw);

%Vectorize time and data so all replicates are used, no need to do this to
%use lsqnonlin
TCAF_10t_all = repmat(TCAF_10t_days, nReps, 1); %should be 45 x 1 
TCAF_10y_all = TCAF_10y_raw(:); %should be 45 x 1 
TCAF_10N = numel(TCAF_10y_all); %total observations so 45

%Define objective function for both models, need to define a residual
%function to be able to use lsqnonlin
logistic_res = @(p, t_dummy) logistic_ode_model(p, TCAF_10t_days, TCAF_10y_raw)-TCAF_10y_all;

genlog_res   = @(p, t_dummy) genlog_ode_model(p, TCAF_10t_days, TCAF_10y_raw)-TCAF_10y_all;

%Define initial guesses for logistic
K0 = 1.4e6; %max(y_all); but it could be higher than the last day counts
r0 = 0.5; %pretty small

p0_log = [K0, r0];
lb_log = [0,  0];     % K >= 0, r >= 0
ub_log = [Inf, Inf]; %could define a differen upper bound

% Create lsqcurvefit problem for logistic
problem_log = createOptimProblem( ...
    'lsqnonlin', ... %eve though we are using a multistart, we still need to define the function to minimize
    'x0',        p0_log, ...
    'objective', logistic_res, ...
    'lb',        lb_log, ...
    'ub',        ub_log);

ms = MultiStart('UseParallel', true,'Display', 'iter', 'StartPointsToRun','bounds-ineqs'); % set false if no Parallel Toolbox
nStarts = 200; %could do more starts if we wanted to 

[p_log, resnorm_log] = run(ms, problem_log, nStarts); %solve the problem with multistart 
fprintf('Fitted logistic parameters for TCAF (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_log(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_log(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_log);
% Fitted logistic parameters for TCAF (10 days)
% Fitted growth rate r = 0.4618
% Fitted carrying capacity k = 163096.7334
% Sum of squared residuals = 8005945209.29

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_log_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(logistic_res, p_log, lb_log, ub_log);

% Compute 95% confidence intervals
ci = nlparci(p_log_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95% Confidence Intervals:\n');
fprintf('r: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
fprintf('K: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% 95r: [0.3979, 0.5257]
% K: [131709.0565, 194491.6834]

%Define initial guesses for generalized logistic
K0  = 1.4e6;         %max(y_all);
r0  = 0.5;
nu0 = 1;            % nu=1 reduces to logistic

p0_gen = [K0, r0, nu0];
lb_gen = [0,  0,   0.1];   % nu > 0
ub_gen = [Inf, Inf, 30];

problem_gen = createOptimProblem( ...
    'lsqnonlin', ...
    'x0',        p0_gen, ...
    'objective', genlog_res, ...
    'lb',        lb_gen, ...
    'ub',        ub_gen);

[p_gen, resnorm_gen] = run(ms, problem_gen, nStarts);
fprintf('Fitted generalized logistic parameters for TCAF (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_gen(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_gen(1));
fprintf('Fitted n term n = %.4f\n', p_gen(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_gen);
% 197 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted generalized logistic parameters for TCAF (10 days)
% Fitted growth rate r = 0.4482
% Fitted carrying capacity k = 344144.3826
% Fitted n term n = 3.2848
% Sum of squared residuals = 9871504074.81

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_gen_final, resnorm_gen_final, residual_gen_final, exitflag_gen, output_gen, lambda_gen, J_gen] = ...
lsqnonlin(genlog_res, p_gen, lb_gen, ub_gen);

% Compute 95% confidence intervals
ci_gen = nlparci(p_gen_final, residual_gen_final, 'jacobian', J_gen);

% Print confidence intervals
fprintf('95% Confidence Intervals:\n');
fprintf('r: [%.4f, %.4f]\n', ci_gen(2,1), ci_gen(2,2));
fprintf('K: [%.4f, %.4f]\n', ci_gen(1,1), ci_gen(1,2));
fprintf('nu: [%.4f, %.4f]\n', ci_gen(3,1), ci_gen(3,2));
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% <stopping criteria details>
% 95r: [0.4252, 0.4711] (0.4711-0.4252)/0.4482 = 0.1024
% K: [316061.3452, 372238.2150] (372238.2150-316061.3452)/344144.3826 =
% 0.1632
% nu: [0.9562, 5.6140] (5.6140-0.9562)/3.2848 = 1.418

%%Plot BIC for the models 
N = numel(TCAF_10y_all);  
k_log = length(p_log);
k_gen = length(p_gen);
SSE_log = resnorm_log; %here the resnorm is already squared
SSE_gen = resnorm_gen;

% BIC formula: N*ln(SSE/N) + k*ln(N)
BIC_log = N * log(SSE_log / N) + k_log * log(N);
BIC_gen = N * log(SSE_gen / N) + k_gen * log(N);

% Display results
fprintf('Logistic:     BIC = %.3f\n', BIC_log); %Logistic:     BIC = 644.122
fprintf('Gen Logistic: BIC = %.3f\n', BIC_gen);

if BIC_log < BIC_gen
    fprintf('Logistic preferred (lower BIC)\n');
else
    fprintf('Gen Logistic preferred (lower BIC)\n');
end

% Logistic:     BIC = 661.814
% Gen Logistic: BIC = 654.531
% Gen Logistic preferred (lower BIC)

%%Plot the curves
% Dense time grid for plotting continuous curves
t_fine = linspace(0, 14, 200)';

% Predict for each replicate separately at t_fine, using fitted params
% Example: use the same ODE solvers but with t_fine instead of 1–14
% (for simplicity, here we interpolate from 1–14)

y_log_fit_full = logistic_ode_model(p_log, TCAF_10t_days, TCAF_10y_raw);
y_gen_fit_full = genlog_ode_model(p_gen, TCAF_10t_days, TCAF_10y_raw);

% Reshape back to 14×3 for plotting
y_log_fit = reshape(y_log_fit_full, [nDays, nReps]);
y_gen_fit = reshape(y_gen_fit_full, [nDays, nReps]);

figure;
hold on;
% Plot data
plot(TCAF_10t_days, TCAF_10y_raw, 'ko', 'MarkerFaceColor',[0.8 0.8 0.8]);

% Plot fitted curves per replicate (logistic)
for j = 1:nReps
    plot(TCAF_10t_days, y_log_fit(:,j), 'b-', 'LineWidth', 1.5);
end

% Plot generalized logistic as dashed red (optional)
for j = 1:nReps
    plot(TCAF_10t_days, y_gen_fit(:,j), 'r--', 'LineWidth', 1.5);
end

xlabel('Time (days)');
ylabel('Cell count (or density)');
legend('Data','Logistic fit','Gen. logistic fit','Location','best');
grid on;
hold off;

saveas(gcf, fullfile(outdir,'TCAF_untreated_logandgenlog_10days.png'));

writematrix(y_log_fit, fullfile(outdir,'TCAF_untreated_log_10days.csv'));
writematrix(y_gen_fit, fullfile(outdir,'TCAF_untreated_genlog_10days.csv'));