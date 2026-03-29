function H = hysteresis_path_bump_polar_cycle()
% hysteresis_path_bump_polar_cycle
% ------------------------------------------------------------
% Custom path hysteresis using chi*.txt recursively.
%
% Paths:
%   path1 : doping start -> end,          polar = 0
%   path2 : doping end   -> doping_bump,  polar = 0
%   path3 : polar  0     -> full,         doping = doping_bump
%   path4 : doping bump  -> start,        polar = full
%   path5 : doping start -> end,          polar = full
%
% Features:
%   - same shifted-chi logic as your working scripts
%   - selectable paths to plot
%   - plots BOTH electronic order psi1 and lattice order X=psi2
%   - plots path map in (doping, polar)
%
% Needs existing functions:
%   - make_realchi_params
%   - make_realchi_coeff
%   - free_energy
%   - minimize_free_energy
% ------------------------------------------------------------

% =========================================================
% 0) params / coeff
% =========================================================
par  = make_realchi_params(true);
coef = make_realchi_coeff(par);

if ~isfield(par, 'plot_paths') || isempty(par.plot_paths)
    par.plot_paths = ["1","2","3","4","5"];
end
if ~isfield(par, 'show_path_map') || isempty(par.show_path_map)
    par.show_path_map = true;
end
if ~isfield(par, 'show_abs_order') || isempty(par.show_abs_order)
    par.show_abs_order = true;
end
if ~isfield(par, 'Npol') || isempty(par.Npol)
    par.Npol = 30;
end

default_root = "/Users/haoranyan/rg_master/data/";
if ~isfolder(default_root)
    default_root = string(pwd);
end
root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
if isequal(root, 0), error('User cancelled.'); end
root = string(root);

fprintf("Root: %s\n", root);

iq_pick = round(par.iq_pick);
jq_pick = round(par.jq_pick);
fprintf("Pick (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);

T_target = par.T_target;
T_tol    = par.T_tol;
fprintf("T_target=%.6g K, T_tol=%.3g K\n", T_target, T_tol);

pmev = double(par.artificial_polar);
dmu_full = 0.5 * pmev * 1e-3;   % meV -> eV, then /2
fprintf("artificial_polar=%.6g meV => dmu_full=%.6g eV\n", pmev, dmu_full);

plot_paths = string(par.plot_paths);

% =========================================================
% 1) collect baseline curve at this T: (mu_i, doping_i, chi_i)
% =========================================================
B = collect_mu_dop_chi_curve_(root, T_target, T_tol, iq_pick, jq_pick);

mu0  = B.mu(:);
dop0 = B.dop(:);
chi0 = B.chi(:);

if isempty(mu0)
    error('No valid chi data found for T_target=%.6g within tol=%.3g at (iq,jq)=(%d,%d).', ...
        T_target, T_tol, iq_pick, jq_pick);
end

fprintf("[baseline] points=%d | mu=[%.6g, %.6g] | dop=[%.6g, %.6g]\n", ...
    numel(mu0), min(mu0), max(mu0), min(dop0), max(dop0));

if numel(mu0) < 5
    error("Too few points at this T. Check T_tol / header parse / (iq,jq).");
end

[mu0, chi0, dop0] = collapse_same_mu_(mu0, chi0, dop0);

Fchi0 = griddedInterpolant(mu0, chi0, "linear", "none");
Fdop0 = griddedInterpolant(mu0, dop0, "linear", "none");

% =========================================================
% 2) define chi(dop, polar=0/full)
% =========================================================
% zero-polar branch
chi_zero_mu = chi0;
ok0 = isfinite(chi_zero_mu);
mu_zero  = mu0(ok0);
chi_zero = chi_zero_mu(ok0);
dop_zero = Fdop0(mu_zero);

ok02 = isfinite(mu_zero) & isfinite(chi_zero) & isfinite(dop_zero);
mu_zero  = mu_zero(ok02);
chi_zero = chi_zero(ok02);
dop_zero = dop_zero(ok02);

[dop_zero, ord0] = sort(dop_zero);
chi_zero = chi_zero(ord0);
mu_zero  = mu_zero(ord0);

[dop_zero, chi_zero, mu_zero] = collapse_same_dop_(dop_zero, chi_zero, mu_zero); %#ok<ASGLU>
chi_of_dop_zero = @(d) interp1(dop_zero, chi_zero, d, 'linear', 'extrap');

% full-polar branch
chi_full_mu = 0.5 * (Fchi0(mu0 + dmu_full) + Fchi0(mu0 - dmu_full));
okf = isfinite(chi_full_mu);
mu_full  = mu0(okf);
chi_full = chi_full_mu(okf);
dop_full = Fdop0(mu_full);

okf2 = isfinite(mu_full) & isfinite(chi_full) & isfinite(dop_full);
mu_full  = mu_full(okf2);
chi_full = chi_full(okf2);
dop_full = dop_full(okf2);

[dop_full, ordf] = sort(dop_full);
chi_full = chi_full(ordf);
mu_full  = mu_full(ordf);

[dop_full, chi_full, mu_full] = collapse_same_dop_(dop_full, chi_full, mu_full); %#ok<ASGLU>
chi_of_dop_full = @(d) interp1(dop_full, chi_full, d, 'linear', 'extrap');

if numel(dop_zero) < 5
    error("Too few zero-polar points after processing.");
end
if numel(dop_full) < 5
    error("Too few full-polar points after shift-average. (dmu too large or mu range too narrow)");
end

% =========================================================
% 3) build base scan window and custom bump path
% =========================================================
[dop_fwd, ~] = build_hyst_doping_grid_(dop_zero, par.hyst);
Nbase = numel(dop_fwd);

dop_start = dop_fwd(1);
dop_end   = dop_fwd(end);

if isfield(par.hyst, 'doping_bump') && isfinite(par.hyst.doping_bump)
    doping_bump = double(par.hyst.doping_bump);
else
    doping_bump = 0.5 * (dop_start + dop_end);
    fprintf('[info] par.hyst.doping_bump missing, use midpoint %.6g\n', doping_bump);
end

% clamp bump into [start,end]
dlo = min(dop_start, dop_end);
dhi = max(dop_start, dop_end);
doping_bump = min(max(doping_bump, dlo), dhi);

fprintf("[custom path] start=%.6g, end=%.6g, bump=%.6g\n", dop_start, dop_end, doping_bump);

dd = median(abs(diff(dop_fwd)));
if ~isfinite(dd) || dd <= 0
    dd = abs(dop_end - dop_start) / max(Nbase-1,1);
end

% path1: start -> end, zero field
dop_1 = dop_fwd(:);

% path2: end -> bump, zero field
dop_2 = build_segment_(dop_end, doping_bump, dd);

% path3: polar ramp at bump
Npol = round(par.Npol);
pol_3 = linspace(0, dmu_full, Npol).';
dop_3 = doping_bump * ones(Npol,1);

% path4: bump -> start, full field
dop_4 = build_segment_(doping_bump, dop_start, dd);

% path5: start -> end, full field
dop_5 = dop_fwd(:);

fprintf("[path sizes] N1=%d, N2=%d, N3=%d, N4=%d, N5=%d\n", ...
    numel(dop_1), numel(dop_2), numel(dop_3), numel(dop_4), numel(dop_5));

scan = struct();
scan.Npsi1 = par.min.Npsi1;
scan.Npsi2 = par.min.Npsi2;
scan.use_gaussian_smooth = par.min.use_gaussian_smooth;
scan.smooth_sigma = par.min.smooth_sigma;
scan.keep_topK = par.min.keep_topK;

opt = struct();
opt.seed_jitter  = par.min.seed_jitter;
opt.seed_repeats = par.min.seed_repeats;
opt.max_iter     = par.min.max_iter;
opt.tol_fun      = par.min.tol_fun;
opt.tol_x        = par.min.tol_x;
opt.cluster_tol  = par.min.cluster_tol;
opt.fd_h         = par.min.fd_h;
opt.min_eig_eps  = par.min.min_eig_eps;
opt.force_origin_and_nearby = par.min.force_origin_and_nearby;
opt.nearby_delta = par.min.nearby_delta;

traj = struct([]);

% =========================================================
% path1: start -> end, polar = 0
% =========================================================
prev = [par.hyst.psi1_0, par.hyst.psi2_0];

psi_1 = nan(numel(dop_1),2);
pol_1 = zeros(numel(dop_1),1);

for i = 1:numel(dop_1)
    dop = dop_1(i);
    chi_used = chi_of_dop_zero(dop);
    C = coef.eval(T_target, dop, chi_used);

    F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
    mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);

    if isempty(mins)
        [p1,p2,~] = minimize_free_energy(C.a1,C.b1,C.a2,C.b2,C.c2,par.lambda, prev(1),prev(2), false);
        pick = [p1,p2];
    else
        pick = pick_closest_min_(mins, prev);
    end

    psi_1(i,:) = pick;
    prev = pick;
end

traj(1).name  = "1";
traj(1).dop   = dop_1;
traj(1).psi   = psi_1;
traj(1).polar = pol_1;

% =========================================================
% path2: end -> bump, polar = 0
% =========================================================
psi_2 = nan(numel(dop_2),2);
pol_2 = zeros(numel(dop_2),1);

for i = 1:numel(dop_2)
    dop = dop_2(i);
    chi_used = chi_of_dop_zero(dop);
    C = coef.eval(T_target, dop, chi_used);

    F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
    mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);

    if isempty(mins)
        [p1,p2,~] = minimize_free_energy(C.a1,C.b1,C.a2,C.b2,C.c2,par.lambda, prev(1),prev(2), false);
        pick = [p1,p2];
    else
        pick = pick_closest_min_(mins, prev);
    end

    psi_2(i,:) = pick;
    prev = pick;
end

traj(2).name  = "2";
traj(2).dop   = dop_2;
traj(2).psi   = psi_2;
traj(2).polar = pol_2;

% =========================================================
% path3: ramp polar at fixed bump doping
% =========================================================
psi_3 = nan(numel(dop_3),2);

for i = 1:numel(dop_3)
    dop = dop_3(i);
    dmu_now = pol_3(i);

    if abs(dmu_now) < 1e-18
        chi_used = chi_of_dop_zero(dop);
    else
        chi_tmp = 0.5 * (Fchi0(mu0 + dmu_now) + Fchi0(mu0 - dmu_now));
        ok_tmp = isfinite(chi_tmp);
        mu_tmp  = mu0(ok_tmp);
        chi_tmp = chi_tmp(ok_tmp);
        dop_tmp = Fdop0(mu_tmp);

        ok_tmp2 = isfinite(mu_tmp) & isfinite(chi_tmp) & isfinite(dop_tmp);
        mu_tmp  = mu_tmp(ok_tmp2);
        chi_tmp = chi_tmp(ok_tmp2);
        dop_tmp = dop_tmp(ok_tmp2);

        [dop_tmp, ord2] = sort(dop_tmp);
        chi_tmp = chi_tmp(ord2);
        mu_tmp  = mu_tmp(ord2);

        [dop_tmp, chi_tmp, mu_tmp] = collapse_same_dop_(dop_tmp, chi_tmp, mu_tmp); %#ok<ASGLU>
        chi_used = interp1(dop_tmp, chi_tmp, dop, 'linear', 'extrap');
    end

    C = coef.eval(T_target, dop, chi_used);

    F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
    mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);

    if isempty(mins)
        [p1,p2,~] = minimize_free_energy(C.a1,C.b1,C.a2,C.b2,C.c2,par.lambda, prev(1),prev(2), false);
        pick = [p1,p2];
    else
        pick = pick_closest_min_(mins, prev);
    end

    psi_3(i,:) = pick;
    prev = pick;
end

traj(3).name  = "3";
traj(3).dop   = dop_3;
traj(3).psi   = psi_3;
traj(3).polar = pol_3;

% =========================================================
% path4: bump -> start, polar = full
% =========================================================
psi_4 = nan(numel(dop_4),2);
pol_4 = dmu_full * ones(numel(dop_4),1);

for i = 1:numel(dop_4)
    dop = dop_4(i);
    chi_used = chi_of_dop_full(dop);
    C = coef.eval(T_target, dop, chi_used);

    F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
    mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);

    if isempty(mins)
        [p1,p2,~] = minimize_free_energy(C.a1,C.b1,C.a2,C.b2,C.c2,par.lambda, prev(1),prev(2), false);
        pick = [p1,p2];
    else
        pick = pick_closest_min_(mins, prev);
    end

    psi_4(i,:) = pick;
    prev = pick;
end

traj(4).name  = "4";
traj(4).dop   = dop_4;
traj(4).psi   = psi_4;
traj(4).polar = pol_4;

% =========================================================
% path5: start -> end, polar = full
% =========================================================
psi_5 = nan(numel(dop_5),2);
pol_5 = dmu_full * ones(numel(dop_5),1);

for i = 1:numel(dop_5)
    dop = dop_5(i);
    chi_used = chi_of_dop_full(dop);
    C = coef.eval(T_target, dop, chi_used);

    F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
    mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);

    if isempty(mins)
        [p1,p2,~] = minimize_free_energy(C.a1,C.b1,C.a2,C.b2,C.c2,par.lambda, prev(1),prev(2), false);
        pick = [p1,p2];
    else
        pick = pick_closest_min_(mins, prev);
    end

    psi_5(i,:) = pick;
    prev = pick;
end

traj(5).name  = "5";
traj(5).dop   = dop_5;
traj(5).psi   = psi_5;
traj(5).polar = pol_5;

% =========================================================
% 4) plot orders vs doping
% =========================================================
if par.show_abs_order
    y11 = abs(psi_1(:,1)); y21 = abs(psi_2(:,1)); y31 = abs(psi_3(:,1)); y41 = abs(psi_4(:,1)); y51 = abs(psi_5(:,1));
    y12 = abs(psi_1(:,2)); y22 = abs(psi_2(:,2)); y32 = abs(psi_3(:,2)); y42 = abs(psi_4(:,2)); y52 = abs(psi_5(:,2));
    ylab1 = '|psi_1|';
    ylab2 = '|X|';
else
    y11 = psi_1(:,1); y21 = psi_2(:,1); y31 = psi_3(:,1); y41 = psi_4(:,1); y51 = psi_5(:,1);
    y12 = psi_1(:,2); y22 = psi_2(:,2); y32 = psi_3(:,2); y42 = psi_4(:,2); y52 = psi_5(:,2);
    ylab1 = '\psi_1';
    ylab2 = 'X';
end

% ----- figure 1: electronic order -----
figure('Color','w'); hold on;
lgd1 = strings(0,1);

h1 = []; h2 = []; h3 = []; h4 = []; h5 = [];

if any(plot_paths == "1")
    h1 = plot(dop_1, y11, '--', 'LineWidth', 2);
    lgd1(end+1) = "path1: fwd";
end
if any(plot_paths == "2")
    h2 = plot(dop_2, y21, '-o', 'LineWidth', 1.8, 'MarkerSize', 4);
    lgd1(end+1) = "path2: to bump";
end
if any(plot_paths == "3")
    h3 = plot(dop_3, y31, ':', 'LineWidth', 2.2);
    lgd1(end+1) = "path3: apply B";
end
if any(plot_paths == "4")
    h4 = plot(dop_4, y41, '-.', 'LineWidth', 2);
    lgd1(end+1) = "path4: to start";
end
if any(plot_paths == "5")
    h5 = plot(dop_5, y51, '-', 'LineWidth', 2.2);
    lgd1(end+1) = "path5: to end";
end

xlabel('doping', 'Interpreter','tex');
ylabel(ylab1, 'Interpreter','tex');
title('Electronic order', 'Interpreter','tex');
grid on; box on;
if ~isempty(lgd1)
    legend(lgd1, 'Location', 'best');
end
if ~isempty(h5)
    uistack(h5, 'top');
end

% ----- figure 2: lattice order X = psi2 -----
figure('Color','w'); hold on;
lgd2 = strings(0,1);

hx1 = []; hx2 = []; hx3 = []; hx4 = []; hx5 = [];

if any(plot_paths == "1")
    hx1 = plot(dop_1, y12, '--', 'LineWidth', 2);
    lgd2(end+1) = "path1: fwd";
end
if any(plot_paths == "2")
    hx2 = plot(dop_2, y22, '-o', 'LineWidth', 1.8, 'MarkerSize', 4);
    lgd2(end+1) = "path2: to bump";
end
if any(plot_paths == "3")
    hx3 = plot(dop_3, y32, ':', 'LineWidth', 2.2);
    lgd2(end+1) = "path3: apply B";
end
if any(plot_paths == "4")
    hx4 = plot(dop_4, y42, '-.', 'LineWidth', 2);
    lgd2(end+1) = "path4: to start";
end
if any(plot_paths == "5")
    hx5 = plot(dop_5, y52, '-', 'LineWidth', 2.2);
    lgd2(end+1) = "path5: to end";
end

xlabel('doping', 'Interpreter','tex');
ylabel(ylab2, 'Interpreter','tex');
title('Lattice order X', 'Interpreter','tex');
grid on; box on;
if ~isempty(lgd2)
    legend(lgd2, 'Location', 'best');
end
if ~isempty(hx5)
    uistack(hx5, 'top');
end

% =========================================================
% 5) plot path map in (doping, polar)
% =========================================================
if par.show_path_map
    figure('Color','w'); hold on;
    lgd3 = strings(0,1);

    if any(plot_paths == "1")
        plot(traj(1).dop, traj(1).polar, '--', 'LineWidth', 2);
        lgd3(end+1) = "path1";
    end
    if any(plot_paths == "2")
        plot(traj(2).dop, traj(2).polar, '-o', 'LineWidth', 1.8, 'MarkerSize', 4);
        lgd3(end+1) = "path2";
    end
    if any(plot_paths == "3")
        plot(traj(3).dop, traj(3).polar, ':', 'LineWidth', 2.2);
        lgd3(end+1) = "path3";
    end
    if any(plot_paths == "4")
        plot(traj(4).dop, traj(4).polar, '-.', 'LineWidth', 2);
        lgd3(end+1) = "path4";
    end
    if any(plot_paths == "5")
        plot(traj(5).dop, traj(5).polar, '-', 'LineWidth', 2.2);
        lgd3(end+1) = "path5";
    end
    
    xlabel('doping', 'Interpreter','tex');
    ylabel('polar shift d\mu (eV)', 'Interpreter','tex');
    title('Scan path in (doping, polar)', 'Interpreter','tex');
    grid on; box on;
    if ~isempty(lgd3)
        legend(lgd3, 'Location', 'best');
    end
end

% =========================================================
% output
% =========================================================
H = struct();
H.traj = traj;
H.baseline.mu  = mu0;
H.baseline.chi = chi0;
H.baseline.dop = dop0;
H.root = root;
H.T_target = T_target;
H.iq_pick = iq_pick;
H.jq_pick = jq_pick;
H.dmu_full = dmu_full;
H.doping_bump = doping_bump;

H.order.path1.psi1 = psi_1(:,1); H.order.path1.X = psi_1(:,2);
H.order.path2.psi1 = psi_2(:,1); H.order.path2.X = psi_2(:,2);
H.order.path3.psi1 = psi_3(:,1); H.order.path3.X = psi_3(:,2);
H.order.path4.psi1 = psi_4(:,1); H.order.path4.X = psi_4(:,2);
H.order.path5.psi1 = psi_5(:,1); H.order.path5.X = psi_5(:,2);

end

% =====================================================================
% Collect (mu, doping, chi) at ONE T and ONE q from many chi*.txt
% =====================================================================
function B = collect_mu_dop_chi_curve_(root, T_target, T_tol, iq_pick, jq_pick)
L = dir(fullfile(root, "**", "chi*.txt"));
L = L(~[L.isdir]);
if isempty(L), error("No chi*.txt under: %s", root); end

Tlo = T_target - T_tol;
Thi = T_target + T_tol;

mu  = nan(numel(L),1);
dop = nan(numel(L),1);
chi = nan(numel(L),1);

n=0; dropT=0; dropH=0; dropTbl=0; dropQ=0;

for k = 1:numel(L)
    fpath  = string(fullfile(L(k).folder, L(k).name));
    folder = string(L(k).folder);

    T_path = parse_T_from_path_(folder);
    H = parse_header_mu_dop_T_(fpath);

    if isfinite(T_path)
        T_use = T_path;
    else
        T_use = H.T;
    end

    if ~isfinite(T_use)
        dropH = dropH + 1;
        continue;
    end

    if ~(T_use >= Tlo && T_use <= Thi)
        dropT = dropT + 1;
        continue;
    end

    mu_path = parse_mu_from_path_(folder);
    mu_use = mu_path;
    if ~isfinite(mu_use)
        mu_use = H.mu;
    end

    if ~isfinite(mu_use) || ~isfinite(H.dop)
        dropH = dropH + 1;
        continue;
    end

    try
        M = readmatrix(fpath, "FileType","text", "CommentStyle","#");
    catch
        dropTbl = dropTbl + 1;
        continue;
    end

    if isempty(M) || size(M,2) < 6
        dropTbl = dropTbl + 1;
        continue;
    end

    Cc = detect_cols_(M);
    id = find(round(M(:,Cc.iq)) == iq_pick & round(M(:,Cc.jq)) == jq_pick, 1);

    if isempty(id)
        dropQ = dropQ + 1;
        continue;
    end

    n = n + 1;
    mu(n)  = mu_use;
    dop(n) = H.dop;
    chi(n) = M(id, Cc.Re);
end

mu  = mu(1:n);
dop = dop(1:n);
chi = chi(1:n);

ok = isfinite(mu) & isfinite(dop) & isfinite(chi);
mu  = mu(ok);
dop = dop(ok);
chi = chi(ok);

fprintf("[collect] files=%d kept=%d | dropT=%d dropH=%d dropTbl=%d dropQ=%d\n", ...
    numel(L), numel(mu), dropT, dropH, dropTbl, dropQ);

B = struct('mu',mu,'dop',dop,'chi',chi);
end

function H = parse_header_mu_dop_T_(fpath)
H = struct('T',NaN,'mu',NaN,'dop',NaN);
fid = fopen(fpath,'r');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

for t = 1:3000
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    s0 = strtrim(string(ln));
    if ~startsWith(s0,"#"), break; end
    s = strtrim(erase(s0,"#"));

    if ~isfinite(H.T)
        tok = regexp(s, "(?:^|\s)T(?:_K)?\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.T = str2double(tok{1}); end
    end

    if ~isfinite(H.dop)
        tok = regexp(s, "(?:^|\s)(?:doping|dop)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.dop = str2double(tok{1}); end
    end

    if ~isfinite(H.mu)
        tok = regexp(s, "(?:^|\s)(?:mu|EF|E_F)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.mu = str2double(tok{1}); end
    end
end
end

function T = parse_T_from_path_(folder)
folder = replace(string(folder), "\", "/");
num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

tok = regexp(folder, "/T\s*=?\s*" + num + "(?:/|$)", "tokens","once","ignorecase");
if ~isempty(tok)
    T = str2double(tok{1});
else
    tok = regexp(folder, "/temp\s*=?\s*" + num + "(?:/|$)", "tokens","once","ignorecase");
    if ~isempty(tok)
        T = str2double(tok{1});
    else
        T = NaN;
    end
end
end

function mu = parse_mu_from_path_(folder)
folder = replace(string(folder), "\", "/");
num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

tok = regexp(folder, "/mu\s*=?\s*" + num + "(?:/|$)", "tokens","once","ignorecase");
if ~isempty(tok)
    mu = str2double(tok{1});
else
    mu = NaN;
end
end

function C = detect_cols_(M)
if size(M,2) >= 8
    C = struct('iq',2,'jq',3,'Re',6);
else
    C = struct('iq',1,'jq',2,'Re',5);
end
end

function [mu2, chi2, dop2] = collapse_same_mu_(mu, chi, dop)
mu  = double(mu(:));
chi = double(chi(:));
dop = double(dop(:));

ok = isfinite(mu) & isfinite(chi) & isfinite(dop);
mu  = mu(ok);
chi = chi(ok);
dop = dop(ok);

[mu_u,~,ic] = unique(mu);
chi_m = accumarray(ic, chi, [], @mean);
dop_m = accumarray(ic, dop, [], @mean);

[mu2,ord] = sort(mu_u);
chi2 = chi_m(ord);
dop2 = dop_m(ord);
end

function [dop2, chi2, mu2] = collapse_same_dop_(dop, chi, mu)
dop = double(dop(:));
chi = double(chi(:));
mu  = double(mu(:));

ok = isfinite(dop) & isfinite(chi) & isfinite(mu);
dop = dop(ok);
chi = chi(ok);
mu  = mu(ok);

[dop_u,~,ic] = unique(dop);
chi_m = accumarray(ic, chi, [], @mean);
mu_m  = accumarray(ic, mu,  [], @mean);

[dop2,ord] = sort(dop_u);
chi2 = chi_m(ord);
mu2  = mu_m(ord);
end

function [dop_fwd, dop_bwd] = build_hyst_doping_grid_(dop0, hyst)
dop0 = double(dop0(:));
dop0 = dop0(isfinite(dop0));

if isempty(dop0)
    error("build_hyst_doping_grid_: dop0 is empty.");
end

dmin_all = min(dop0);
dmax_all = max(dop0);

mode = "N";
if isfield(hyst, "grid_mode") && ~isempty(hyst.grid_mode)
    mode = string(hyst.grid_mode);
end

direction = "ascend";
if isfield(hyst, "forward_direction") && ~isempty(hyst.forward_direction)
    direction = lower(string(hyst.forward_direction));
end

if ~(direction == "ascend" || direction == "descend")
    error("hyst.forward_direction must be 'ascend' or 'descend'");
end

scan_start = NaN;
scan_end   = NaN;

if isfield(hyst, "scan_start") && ~isempty(hyst.scan_start)
    scan_start = double(hyst.scan_start);
end
if isfield(hyst, "scan_end") && ~isempty(hyst.scan_end)
    scan_end = double(hyst.scan_end);
end

if ~isfinite(scan_start), scan_start = dmin_all; end
if ~isfinite(scan_end),   scan_end   = dmax_all; end

dlo_req = min(scan_start, scan_end);
dhi_req = max(scan_start, scan_end);

dlo = max(dmin_all, dlo_req);
dhi = min(dmax_all, dhi_req);

if ~(isfinite(dlo) && isfinite(dhi) && dhi > dlo)
    error("Requested scan window [%.6g, %.6g] has no overlap with data range [%.6g, %.6g].", ...
        dlo_req, dhi_req, dmin_all, dmax_all);
end

fprintf("[build grid] data range = [%.6g, %.6g], requested = [%.6g, %.6g], used = [%.6g, %.6g]\n", ...
    dmin_all, dmax_all, dlo_req, dhi_req, dlo, dhi);

if mode == "step"
    step = double(hyst.step);
    if ~isfinite(step) || step <= 0
        error("hyst.step must be > 0 when grid_mode='step'");
    end

    dop_base = (dlo:step:dhi).';
    if isempty(dop_base) || abs(dop_base(end) - dhi) > 1e-12
        dop_base(end+1,1) = dhi; %#ok<AGROW>
    end

elseif mode == "N"
    if isfield(hyst, "N")
        N = double(hyst.N);
    else
        error("hyst.N is required when grid_mode='N'");
    end
    if ~isfinite(N) || N < 2
        error("hyst.N must be >= 2 when grid_mode='N'");
    end
    dop_base = linspace(dlo, dhi, round(N)).';

else
    error("Unsupported hyst.grid_mode. Use 'N' or 'step'.");
end

if direction == "ascend"
    dop_fwd = dop_base;
else
    dop_fwd = flipud(dop_base);
end

dop_bwd = flipud(dop_fwd);
end

function seg = build_segment_(x0, x1, dx)
if abs(x1 - x0) < 1e-14
    seg = x0;
    return;
end

if ~isfinite(dx) || dx <= 0
    seg = [x0; x1];
    return;
end

n = max(2, round(abs(x1 - x0) / dx) + 1);
seg = linspace(x0, x1, n).';
end

function pick = pick_closest_min_(mins, prev)
d = sqrt(sum((mins - prev).^2, 2));
[~,ix] = min(d);
pick = mins(ix,:);
end

function mins = find_local_minima_2d_denseSeeds_(F, psi1_lim, psi2_lim, scan, opt)
psi1 = linspace(psi1_lim(1), psi1_lim(2), scan.Npsi1);
psi2 = linspace(psi2_lim(1), psi2_lim(2), scan.Npsi2);
[P1,P2] = meshgrid(psi1, psi2);

Fgrid = arrayfun(@(u,v) safe_F_(F,u,v,psi1_lim,psi2_lim), P1, P2);

if scan.use_gaussian_smooth
    Fgrid = gaussian_smooth_2d_(Fgrid, scan.smooth_sigma);
end

mask_min = discrete_local_min_mask_(Fgrid);
[rr,cc] = find(mask_min);

if isempty(rr)
    seeds = [0,0];
else
    seeds = [P1(sub2ind(size(P1), rr, cc)), P2(sub2ind(size(P2), rr, cc))];
    E = Fgrid(sub2ind(size(Fgrid), rr, cc));
    [~,ix] = sort(E,'ascend');
    ix = ix(1:min(scan.keep_topK, numel(ix)));
    seeds = seeds(ix,:);
end

if isfield(opt,'force_origin_and_nearby') && opt.force_origin_and_nearby
    d = opt.nearby_delta;
    seeds = [seeds;
             0,0;
             d,0; -d,0; 0,d; 0,-d];
end

obj = @(x) safe_F_(F, x(1), x(2), psi1_lim, psi2_lim);

fopt = optimset('Display','off', ...
    'MaxIter', opt.max_iter, ...
    'TolFun',  opt.tol_fun, ...
    'TolX',    opt.tol_x);

cand = zeros(0,2);
for k = 1:size(seeds,1)
    x0 = seeds(k,:);
    for rep = 1:opt.seed_repeats
        xstart = x0 + opt.seed_jitter * randn(1,2);
        xstart(1) = min(max(xstart(1), psi1_lim(1)), psi1_lim(2));
        xstart(2) = min(max(xstart(2), psi2_lim(1)), psi2_lim(2));
        try
            xhat = fminsearch(obj, xstart, fopt);
        catch
            continue;
        end
        xhat(1) = min(max(xhat(1), psi1_lim(1)), psi1_lim(2));
        xhat(2) = min(max(xhat(2), psi2_lim(1)), psi2_lim(2));
        cand(end+1,:) = xhat; %#ok<AGROW>
    end
end

if isempty(cand)
    mins = zeros(0,2);
    return;
end

uniq = cluster_points_(cand, opt.cluster_tol);

keep = false(size(uniq,1),1);
for i = 1:size(uniq,1)
    x = uniq(i,:);
    Hh = hessian_fd_( @(u,v) F(u,v), x(1), x(2), opt.fd_h );
    ev = eig((Hh+Hh')/2);
    if all(real(ev) > opt.min_eig_eps)
        keep(i) = true;
    end
end

mins = uniq(keep,:);
end

function val = safe_F_(F, psi1, psi2, lim1, lim2)
if psi1 < lim1(1) || psi1 > lim1(2) || psi2 < lim2(1) || psi2 > lim2(2)
    val = 1e30;
    return;
end
val = F(psi1, psi2);
if ~isfinite(val)
    val = 1e30;
end
end

function uniq = cluster_points_(P, tol)
uniq = zeros(0,2);
for i = 1:size(P,1)
    x = P(i,:);
    if isempty(uniq)
        uniq = x;
    else
        d = sqrt(sum((uniq - x).^2, 2));
        if all(d > tol)
            uniq(end+1,:) = x; %#ok<AGROW>
        end
    end
end
end

function Hh = hessian_fd_(f, x, y, h)
f00 = f(x,y);
fxx = (f(x+h,y) - 2*f00 + f(x-h,y))/h^2;
fyy = (f(x,y+h) - 2*f00 + f(x,y-h))/h^2;
fxy = (f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)) / (4*h^2);
Hh = [fxx, fxy; fxy, fyy];
end

function mask = discrete_local_min_mask_(A)
[R,C] = size(A);
mask = false(R,C);
if R < 3 || C < 3, return; end

center = A(2:R-1, 2:C-1);
n1 = A(1:R-2, 2:C-1);
n2 = A(3:R  , 2:C-1);
n3 = A(2:R-1, 1:C-2);
n4 = A(2:R-1, 3:C  );
n5 = A(1:R-2, 1:C-2);
n6 = A(1:R-2, 3:C  );
n7 = A(3:R  , 1:C-2);
n8 = A(3:R  , 3:C  );

m = center <= n1 & center <= n2 & center <= n3 & center <= n4 & ...
    center <= n5 & center <= n6 & center <= n7 & center <= n8;

mask(2:R-1, 2:C-1) = m;
end

function B = gaussian_smooth_2d_(A, sigma)
if sigma <= 0
    B = A;
    return;
end
rad = max(1, ceil(3*sigma));
x = (-rad:rad);
g = exp(-(x.^2)/(2*sigma^2));
g = g / sum(g);
B = conv2(A, g, 'same');
B = conv2(B, g', 'same');
end