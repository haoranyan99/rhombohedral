function coef = make_realchi_coeff(par)
% Return coefficient evaluator using ONLY (T, dop, chi_used) and par.
    Tc = par.Xboundary_Tc;
    dc = par.Xboundary_dc;
    Tm = par.Xboundary_Tmin;
    dm = par.Xboundary_dmin;
    a2 = par.Xboundary_a2_T;

    coef = struct();
    coef.beta = @(T) par.Xboundary_Tslope.*T + par.Xboundary_beta0; 

    coef.b1 = @(T,dop) 1;
    coef.c2 = @(T,dop) 1;
    coef.b2 = @(T) coef.beta(T);
    coef.a2 = @(T,dop) 5.0 / 6.0 * (coef.beta(T).^2) ...
        - a2*(T-Tc)^2 + a2*(Tm-Tc)^2*(dc-dop)/(dc-dm);

    coef.eval = @(T,dop,chi_used) eval_all_(par, coef, T, dop, chi_used);
end

function C = eval_all_(par, coef, T, dop, chi_used)
    C = struct();
    C.a1 = par.chi_scaler * (par.invV + abs(chi_used));
    C.b1 = coef.b1(T,dop);
    C.a2 = coef.a2(T,dop);
    C.b2 = coef.b2(T);
    C.c2 = coef.c2(T,dop);
end
