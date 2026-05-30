function coef = make_realchi_coeff(par)

    coef = struct();
    c2 = 1; % c2 > 0 to ensure F>0 when X is large

    % calculate b2 from theoretical boundary 5/6 * b2^2 - a2 * c2 = 0

    % a2 = alpha * (T - Tc)
    alpha = par.Lat_alpha;
    Tc = par.Lat_Tc;

    coef.a2 = @(T) alpha * (T-Tc);
    coef.T_boundary = @(dop) lattice_boundary_T_(par, dop);
    coef.b2 = @(dop) -sqrt(max(0, 1.2 * c2 * alpha * (coef.T_boundary(dop) - Tc)));
    coef.c2 = c2;

    coef.eval = @(T,dop,chi_used) eval_all_(par, coef, T, dop, chi_used);
end

function Tbd = lattice_boundary_T_(par, dop)

mode = "parabola";
if isfield(par, "Lat_boundary_mode")
    mode = string(par.Lat_boundary_mode);
end

switch lower(mode)
    case "sigmoid_step"
        Tbd = lattice_boundary_T_sigmoid_step_(par, dop);

    otherwise
        dop1 = par.Lat_pt1(1);
        T1 = par.Lat_pt1(2);
        dop2 = par.Lat_pt2(1);
        T2 = par.Lat_pt2(2);
        A = (T1 - T2) / (dop2*dop2 - dop1*dop1);
        C = (T1*dop2*dop2 - T2*dop1*dop1) / (dop2*dop2 - dop1*dop1);
        Tbd = -A .* dop .* dop + C;
end

end

function Tbd = lattice_boundary_T_sigmoid_step_(par, dop)

s = par.Lat_step;
Tmin = s.T_min_K;
Tk = s.T_knee_K;
Tmax = s.T_max_K;
xL = s.dop_left;
xK = s.dop_knee;
xR = s.dop_right;
p = s.slow_power;
w = s.fast_width_K;
n_grid = s.n_grid;

Tgrid = linspace(Tmin, Tmax, n_grid);
xgrid = nan(size(Tgrid));

idx_slow = Tgrid <= Tk;
u = (Tgrid(idx_slow) - Tmin) ./ max(Tk - Tmin, eps);
u = min(max(u, 0), 1);
xgrid(idx_slow) = xL + (xK - xL) .* u.^p;

idx_fast = ~idx_slow;
v = sigmoid01_((Tgrid(idx_fast) - Tk) ./ max(w, eps));
v0 = sigmoid01_(0);
v1 = sigmoid01_((Tmax - Tk) ./ max(w, eps));
v = (v - v0) ./ max(v1 - v0, eps);
v = min(max(v, 0), 1);
xgrid(idx_fast) = xK + (xR - xK) .* v;

[xuniq, ia] = unique(xgrid, "stable");
Tuniq = Tgrid(ia);

Tbd = interp1(xuniq, Tuniq, dop, "linear", "extrap");
Tbd = min(max(Tbd, Tmin), Tmax);

end

function y = sigmoid01_(x)

y = 1 ./ (1 + exp(-x));

end

function C = eval_all_(par, coef, T, dop, chi_used)
    C = struct();
    C.a1 = par.chi_scaler * (par.invV - abs(chi_used));
    C.b1 = 1;
    C.a2 = coef.a2(T);
    C.b2 = coef.b2(dop);
    % C.b2 = -5;
    C.c2 = coef.c2;
end
