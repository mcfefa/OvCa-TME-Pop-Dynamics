function [res, BIC_conditions] = coculture_residuals_3pop_untreated_model1(p, rS, K, alphaS, IC50S, rR, alphaR, IC50R, rF, KF, t_data, S_raw_all, R_raw_all)

    betaS = p(1);
    betaR = p(2);

    res = [];
    BIC_conditions = zeros(3,1);

    for i = 1:3
        S_raw = S_raw_all{i};
        R_raw = R_raw_all{i};
        [~, nReps] = size(S_raw);

        res_condition = [];

        for j = 1:nReps
            S0 = S_raw(1,j);
            R0 = R_raw(1,j);
            F0 = 10000;

            [S_mod, R_mod, ~] = simulate_coculture_3pop_untreated_model1( ...
                t_data, S0, R0, F0, p, rS, K, alphaS, IC50S, rR, alphaR, IC50R, rF, KF);

            resS = S_raw(:,j) - S_mod;
            resR = R_raw(:,j) - R_mod;

            res_condition = [res_condition; resS; resR];
        end

        res = [res; res_condition];

        N_cond = length(res_condition);
        k = length(p);
        SSE_cond = sum(res_condition.^2);
        BIC_conditions(i) = N_cond * log(SSE_cond / N_cond) + k * log(N_cond);
    end
end