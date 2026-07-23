function dydt = coculture_ode_3pop_untreated_model1(t, y, rS, K, alphaS, IC50S, rR, alphaR, IC50R, rF, KF, betaS, betaR)

    S = y(1);
    R = y(2);
    F = y(3);

    dSdt = rS * S * (1 - (S + R)/K) + betaS * F;
    dRdt = rR * R * (1 - (R + S)/K) + betaR * F;
    dFdt = rF * F * (1 - F/KF);

    dydt = [dSdt; dRdt; dFdt];
end