function dydt = coculture_ode_3pop_untreated_model2(t, y, rS, K, alphaS, IC50S, rR, alphaR, IC50R)

    S = y(1);
    R = y(2);

    dSdt = rS * S * (1 - (S + R)/K);
    dRdt = rR * R * (1 - (R + S)/K);

    dydt = [dSdt; dRdt];
end