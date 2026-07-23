function dydt = coculture_ode_adaptive(t, y, p, D_t)
    S = y(1);
    R = y(2);
    F = y(3);

    kill_S = (D_t) ./ (p.IC50S + D_t);
    kill_R = (D_t) ./ (p.IC50R + D_t);

    alpha_tS = p.alphaS * t;
    alpha_tR = p.alphaR * t;

    dSdt = p.rS * S * (1 - (S + R - p.betaS * F)/p.K) - alpha_tS * kill_S * S;
    dRdt = p.rR * R * (1 - (R + S - p.betaR * F)/p.K) - alpha_tR * kill_R * R;
    dFdt = p.rF * F * (1 - F/p.KF);

    dydt = [dSdt; dRdt; dFdt];
end