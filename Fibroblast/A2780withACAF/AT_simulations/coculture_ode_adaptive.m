function dydt = coculture_ode_adaptive(t, y, p, D_t)
    S = y(1);
    R = y(2);
    F = y(3);

    kill_S = (D_t.^p.nS) ./ (p.IC50S.^p.nS + D_t.^p.nS);
    kill_R = (D_t.^p.nR) ./ (p.IC50R.^p.nR + D_t.^p.nR);

    alpha_tS = p.alphaS * t;
    alpha_tR = p.alphaR * t;

    dSdt = p.rS * S * (1 - (S + R)/p.K) - alpha_tS * kill_S * S + p.betaS * F;
    dRdt = p.rR * R * (1 - (R + S)/p.K) - alpha_tR * kill_R * R - p.betaR * F;
    dFdt = p.rF * F * (1 - F/p.KF);

    dydt = [dSdt; dRdt; dFdt];
end