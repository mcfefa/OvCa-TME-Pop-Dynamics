function dydt = coculture_ode_sens2_m13(t, y, rS, K, alphaS, nS, IC50S, rR, alphaR, nR, IC50R, ce, betaS, betaR, rF, KF)
    S = y(1); R = y(2); F = y(3);
    D_t = constant_D_schedule(t, ce);  % always = Ce
    %hill kill terms, different for S (alpha is time dependent) and R
    kill_S = (D_t.^nS) ./ (IC50S.^nS + D_t.^nS);
    alpha_tS = alphaS * t;
    kill_R = (D_t.^nR) ./ (IC50R.^nR + D_t.^nR);
    alpha_tR = alphaR * t; 

    dSdt = rS * S * (1 - (S + R)/K) - alpha_tS * kill_S * S + betaS * F;
    dRdt = rR * R * (1 - (R + S)/K) - alpha_tR * kill_R * R - betaR * F;
    dFdt = rF * F * (1 - F/KF);
    dydt = [dSdt; dRdt; dFdt];
end