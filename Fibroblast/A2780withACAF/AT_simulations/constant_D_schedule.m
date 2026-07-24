function D_t = constant_D_schedule(t_day, ce)
    % Constant drug concentration at ALL time points
    D_t = ce * ones(size(t_day));
end