function coef = make_realchi_coeff(par)

    coef = struct();
    c2 = 1; % c2 > 0 to ensure F>0 when X is large

    % new def: T = -A * dop^2 + C for lattice phase boundary
    dop1 = par.Lat_pt1(1);
    T1 = par.Lat_pt1(2);
    dop2 = par.Lat_pt2(1);
    T2 = par.Lat_pt2(2);
    A = (T1 - T2) / (dop2*dop2 - dop1*dop1);
    C = (T1*dop2*dop2 - T2*dop1*dop1) / (dop2*dop2 - dop1*dop1);
    

    % calculate b2 from theoretical boundary 5/6 * b2^2 - a2 * c2 = 0

    % a2 = alpha * (T - Tc)
    alpha = par.Lat_alpha;
    Tc = par.Lat_Tc;

    coef.a2 = @(T) alpha * (T-Tc);
    coef.b2 = @(dop) sqrt(max(0, 1.2 * c2 * alpha * (-A * dop * dop + C - Tc)));
    coef.c2 = c2;

    coef.eval = @(T,dop,chi_used) eval_all_(par, coef, T, dop, chi_used);
end

function C = eval_all_(par, coef, T, dop, chi_used)
    C = struct();
    C.a1 = par.chi_scaler * (par.invV + abs(chi_used));
    C.b1 = 1;
    C.a2 = coef.a2(T);
    C.b2 = coef.b2(dop);
    % C.b2 = -5;
    C.c2 = coef.c2;
end
