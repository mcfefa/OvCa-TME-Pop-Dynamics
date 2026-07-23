function dydt = coculture_ode_3pop_untreated_model3(t, y, rS, K, alphaS, nS, IC50S, rR, alphaR, nR, IC50R, rF, KF, betaS, betaR)

    S = y(1);
    R = y(2);
    F = y(3);

    dSdt = rS * S * (1 - (S + R - betaS * F)/K);
    dRdt = rR * R * (1 - (R + S + betaR * F)/K);
    dFdt = rF * F * (1 - F/KF);

    dydt = [dSdt; dRdt; dFdt];
end