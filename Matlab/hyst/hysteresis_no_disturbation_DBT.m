function H = hysteresis_no_disturbation_DBT()

this_dir = fileparts(mfilename("fullpath"));
addpath(this_dir, "-begin");
output_dir = fullfile(this_dir, "hyst_data");
if ~isfolder(output_dir)
    mkdir(output_dir);
end
run_stamp = string(datetime("now", "Format", "yyyy_MM_dd_HHmmss"));

par  = make_realchi_params(true);
coef = make_realchi_coeff(par);

% =========================
% USER SETTINGS
% =========================
% Current chi tree:
%   root/D-0.010/B0T/T6.400/mu0.800000/chi_valley_..._q5_3_...txt
root = "E:/rg_master/rhombohedral/data";
D_folder = "D-0.010";

% [] means auto-detect all available B/T folders under root/D_folder.
B_list = [0];
T_list = [6.4];

% "average" averages all chi_valley_*_spin_* files at the same mu.
% Or use one branch, e.g. "valley_plus_spin_up".
chi_branch = "average";
doping_scale = 2.0;
broadening_meV = 0.1;  % [] means use all broadenings found in filenames.

use_parallel = true;   % <<< 是否并行
par.hyst.force_fwd_psi_zero = true;
par.hyst.force_fwd_psi_zero_window = [-1.55, -1.45];
par.io.doping_scale = doping_scale;
par.io.broadening_meV = broadening_meV;

plot_a1  = false;
plot_psi = false;
plot_X   = false;

plot_a1  = true;
plot_psi = true;
plot_X   = true;

if ~isfolder(root)
    default_root = "E:/rg_master/rhombohedral/data/";
    if ~isfolder(default_root)
        default_root = "/Users/haoranyan/rg_master/rhombohedral/data/";
    end
    if ~isfolder(default_root)
        default_root = pwd;
    end

    root = uigetdir(default_root, 'Select chi root folder');
    if isequal(root,0), error('User cancelled.'); end
    root = string(root);
end

fprintf("Root: %s\n", root);
fprintf("D folder: %s\n", D_folder);
fprintf("Output dir: %s\n", output_dir);
fprintf("Run stamp: %s\n", run_stamp);
fprintf("chi branch: %s\n", chi_branch);
if isempty(broadening_meV)
    fprintf("broadening: all\n");
else
    fprintf("broadening: %.12g meV\n", broadening_meV);
end

iq_pick = round(par.iq_pick);
jq_pick = round(par.jq_pick);

fprintf("Pick (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);

BTfolders = detect_BT_folders_(root, D_folder);
if isempty(BTfolders)
    error("No D/B/T folders found under root: %s", root);
end

if isempty(B_list)
    B_list = unique([BTfolders.B].', 'stable');
else
    B_list = B_list(:);
end

if isempty(T_list)
    T_list = unique([BTfolders.T].', 'stable');
else
    T_list = T_list(:);
end

nB = numel(B_list);
nT = numel(T_list);

Res = cell(nB, nT);

fprintf("B list: %s T\n", strjoin(compose("%.6g", B_list), ", "));
fprintf("T list: %s K\n", strjoin(compose("%.6g", T_list), ", "));

% =========================
% Build task list
% =========================
tasks = [];
for ib = 1:nB
    for it = 1:nT
        tasks(end+1,:) = [ib, it]; %#ok<AGROW>
    end
end

nTask = size(tasks,1);
ResTask = cell(nTask,1);

% =========================
% Parallel or serial run
% =========================
if use_parallel && nTask > 1

    pool = gcp('nocreate');
    if isempty(pool)
        parpool;
    end

    parfor kk = 1:nTask

        ib = tasks(kk,1);
        it = tasks(kk,2);

        B_target = B_list(ib);
        T_target = T_list(it);
        dmu = 0.0;

        T_folder = find_BT_folder_(BTfolders, B_target, T_target, par.T_tol);

        if strlength(T_folder) == 0
            warning("Skip B=%.6g T, T=%.6g K: no matched folder.", B_target, T_target);
            ResTask{kk} = [];
            continue;
        end

        fprintf("\n[parfor] B = %.6g T, T = %.6g K\n", ...
            B_target, T_target);

        try
            R = run_one_T_original_minima_(T_folder, T_target, ...
                iq_pick, jq_pick, dmu, par, coef, chi_branch, doping_scale, broadening_meV);
            R.B_T = B_target;
            R.B_folder = string(fileparts(char(T_folder)));
            R.D_folder = D_folder;
            R.chi_branch = chi_branch;
        catch ME
            warning("Skip B=%.6g T, T=%.6g K: %s", ...
                B_target, T_target, ME.message);
            R = [];
        end

        ResTask{kk} = R;
    end

    for kk = 1:nTask
        ib = tasks(kk,1);
        it = tasks(kk,2);
        Res{ib,it} = ResTask{kk};
    end

else

    for kk = 1:nTask

        ib = tasks(kk,1);
        it = tasks(kk,2);

        B_target = B_list(ib);
        T_target = T_list(it);
        dmu = 0.0;

        T_folder = find_BT_folder_(BTfolders, B_target, T_target, par.T_tol);

        if strlength(T_folder) == 0
            warning("Skip B=%.6g T, T=%.6g K: no matched folder.", B_target, T_target);
            continue;
        end

        fprintf("\n================ B = %.6g T, T = %.6g K ================\n", ...
            B_target, T_target);
        fprintf("Use T folder: %s\n", T_folder);
        fprintf("---- real magnetic chi, dmu = %.6g eV ----\n", dmu);

        try
            R = run_one_T_original_minima_(T_folder, T_target, ...
                iq_pick, jq_pick, dmu, par, coef, chi_branch, doping_scale, broadening_meV);
            R.B_T = B_target;
            R.B_folder = string(fileparts(char(T_folder)));
            R.D_folder = D_folder;
            R.chi_branch = chi_branch;
        catch ME
            warning("Skip B=%.6g T, T=%.6g K: %s", ...
                B_target, T_target, ME.message);
            continue;
        end

        Res{ib, it} = R;
    end
end

% =========================
% Serial save and plot
% =========================
for ib = 1:nB
    for it = 1:nT
            R = Res{ib,it};

            if isempty(R)
                continue;
            end

            write_hyst_data_BT(R, ib, it, nB, nT, output_dir, run_stamp);

            plot_one_T_(R, par, plot_a1, plot_psi, plot_X);
    end
end

H = struct();
H.root = root;
H.D_folder = D_folder;
H.B_list = B_list;
H.T_list = T_list;
H.chi_branch = chi_branch;
H.doping_scale = doping_scale;
H.broadening_meV = broadening_meV;
H.results = Res;
H.iq_pick = iq_pick;
H.jq_pick = jq_pick;
H.par = par;
H.output_dir = output_dir;
H.run_stamp = run_stamp;

end

function R = run_one_T_original_minima_(T_folder, T_target, iq_pick, jq_pick, dmu, par, coef, chi_branch, doping_scale, broadening_meV)

B = collect_mu_dop_chi_curve_(T_folder, iq_pick, jq_pick, chi_branch, broadening_meV);

mu0  = B.mu(:);
dop0 = B.dop(:);
chi0 = B.chi(:);
dop_raw_from_file = dop0;
dop0 = doping_scale * dop0;

fprintf("[baseline] points=%d | mu=[%.6g, %.6g] | dop_raw=[%.6g, %.6g] | dop_scaled=[%.6g, %.6g]\n", ...
    numel(mu0), min(mu0), max(mu0), ...
    min(dop_raw_from_file), max(dop_raw_from_file), min(dop0), max(dop0));

if numel(mu0) < 5
    error("Too few points in T folder. Check header parse / (iq,jq).");
end

if abs(dmu) > 1e-18
    warning("dmu is ignored in the real-B DBT hysteresis path. Use the artificial-polar script for mu-shift averaging.");
end

ok = isfinite(mu0) & isfinite(dop0) & isfinite(chi0);
mu_use  = mu0(ok);
dop_use = dop0(ok);
chi_use = chi0(ok);

if numel(mu_use) < 5
    error("Too few finite chi points after filtering.");
end

% Match hysteresis_nodisturbation_realPolar: build chi directly on the
% doping axis, then average only duplicate doping points.
[dop_use, ord] = sort(dop_use);
chi_use = chi_use(ord);
mu_use  = mu_use(ord);

[dop_use, chi_use, mu_use] = collapse_same_dop_(dop_use, chi_use, mu_use);

chi_unsmoothed = chi_use;
[chi_use, chi_smooth_info] = smooth_chi_input_curve_(dop_use, chi_use, par.smooth);

chi_of_dop = @(d) interp1(dop_use, chi_use, d, 'linear', 'extrap');

[dop_fwd, dop_bwd] = build_hyst_doping_grid_(dop_use, par.hyst);
N = numel(dop_fwd);

fprintf("[scan window] forward: [%.6g, %.6g], backward: [%.6g, %.6g], N=%d\n", ...
    dop_fwd(1), dop_fwd(end), dop_bwd(1), dop_bwd(end), N);

scan = make_scan_options_(par);
opt = make_min_options_(par);

psi_f = nan(N,2);
a1_f  = nan(N,1);
chi_f = nan(N,1);

prev = [par.hyst.psi1_0, par.hyst.psi2_0];

for i = 1:N
    dop = dop_fwd(i);
    chi_used = chi_of_dop(dop);
    C = coef.eval(T_target, dop, chi_used);

    pick = choose_local_minimum_(C, par, scan, opt, prev, true);

    if isfield(par.hyst, "force_fwd_psi_zero") && par.hyst.force_fwd_psi_zero
        win = par.hyst.force_fwd_psi_zero_window;
        if dop >= min(win) && dop <= max(win)
            pick(1) = 0;
        end
    end

    psi_f(i,:) = pick;
    prev = pick;

    a1_f(i)  = C.a1;
    chi_f(i) = chi_used;
end

psi_b = nan(N,2);
a1_b  = nan(N,1);
chi_b = nan(N,1);

prev = psi_f(end,:);

for i = 1:N
    dop = dop_bwd(i);
    chi_used = chi_of_dop(dop);
    C = coef.eval(T_target, dop, chi_used);

    pick = choose_local_minimum_(C, par, scan, opt, prev, false);

    psi_b(i,:) = pick;
    prev = pick;

    a1_b(i)  = C.a1;
    chi_b(i) = chi_used;
end


% ============================================================
% Post-check for backward path:
% Run backward minimization first, then enforce:
% if current a1 > 0 and previous point is already trivial,
% force current point to remain at (psi,X)=(0,0).
% This is part of the no-disturbance hysteresis definition.
% ============================================================
trivial_tol = 1e-3;

for i = 2:N

    prev_is_trivial = ...
        abs(psi_b(i-1,1)) < trivial_tol && ...
        abs(psi_b(i-1,2)) < trivial_tol;

    if a1_b(i) > 0 && prev_is_trivial
        psi_b(i,:) = [0, 0];
    end

end

psi_f_abs = abs(psi_f);
psi_b_abs = abs(psi_b);

psi_f_abs_plot = psi_f_abs;
psi_b_abs_plot = psi_b_abs;
gamma_used = NaN;

R = struct();
R.T = T_target;
R.T_folder = T_folder;
R.iq_pick = iq_pick;
R.jq_pick = jq_pick;
R.dmu_eV = dmu;
R.doping_scale = doping_scale;
R.broadening_meV = broadening_meV;

R.mu_raw  = mu0;
R.dop_raw_from_file = dop_raw_from_file;
R.dop_raw = dop0;
R.chi_raw = chi0;

R.mu_used  = mu_use;
R.dop_used = dop_use;
R.chi_used = chi_use;
R.chi_unsmoothed = chi_unsmoothed;
R.chi_smooth_info = chi_smooth_info;

R.dop_fwd = dop_fwd;
R.dop_bwd = dop_bwd;

R.psi_fwd = psi_f;
R.psi_bwd = psi_b;

R.psi_f_abs_plot = psi_f_abs_plot;
R.psi_b_abs_plot = psi_b_abs_plot;

R.a1_fwd  = a1_f;
R.a1_bwd  = a1_b;

R.chi_fwd = chi_f;
R.chi_bwd = chi_b;

R.gamma_used = gamma_used;
R.par = par;
R.run_params = make_run_params_record_(R, par);

end

function run_params = make_run_params_record_(R, par)

run_params = struct();
run_params.created_by = "hysteresis_no_disturbation_DBT";
run_params.created_at = string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss"));
run_params.T = R.T;
run_params.B_T = getfield_or_(R, "B_T", NaN);
run_params.D_folder = getfield_or_(R, "D_folder", "");
run_params.T_folder = R.T_folder;
run_params.chi_branch = getfield_or_(R, "chi_branch", "");
run_params.iq_pick = R.iq_pick;
run_params.jq_pick = R.jq_pick;
run_params.doping_scale = R.doping_scale;
run_params.broadening_meV = R.broadening_meV;
run_params.chi_smooth_info = R.chi_smooth_info;
run_params.doping_note = "doping in chi headers is multiplied by doping_scale before chi(doping) interpolation and hysteresis.";
run_params.chi_to_a1 = "a1 = par.chi_scaler * (par.invV - abs(chi))";
run_params.lattice_boundary = "T_boundary(dop) is selected by par.Lat_boundary_mode; sigmoid_step inverts the configured x_c(T) boundary.";
run_params.lattice_b2 = "b2(dop) = -sqrt(max(0, 1.2*c2*par.Lat_alpha*(T_boundary(dop)-par.Lat_Tc))) from delta=(5/6)b2^2-a2*c2=0.";
run_params.par = par;

end

function val = getfield_or_(S, name, fallback)

if isfield(S, name)
    val = S.(name);
else
    val = fallback;
end

end

function pick = choose_local_minimum_(C, par, scan, opt, prev, allow_zero_if_a1_positive)

if allow_zero_if_a1_positive && C.a1 > 0
    pick = [0, 0];
    return;
end

F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);

if isempty(mins)
    [p1,p2,~] = minimize_free_energy( ...
        C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda, ...
        prev(1), prev(2), false);
    pick = [p1, p2];
else
    pick = pick_closest_min_(mins, prev);
end

end

function scan = make_scan_options_(par)

scan = struct();
scan.Npsi1 = par.min.Npsi1;
scan.Npsi2 = par.min.Npsi2;
scan.use_gaussian_smooth = par.min.use_gaussian_smooth;
scan.smooth_sigma = par.min.smooth_sigma;
scan.keep_topK = par.min.keep_topK;

end

function opt = make_min_options_(par)

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

end

function [chi_out, info] = smooth_chi_input_curve_(dop, chi_in, smooth)

chi_out = chi_in;
info = struct();
info.enabled = false;
info.method = "none";
info.gamma = NaN;

if ~isfield(smooth, "use_lorentz") || ~smooth.use_lorentz
    return;
end

gamma = NaN;
if isfield(smooth, "lorentz_gamma")
    gamma = smooth.lorentz_gamma;
end

if ~isfinite(gamma) || gamma <= 0
    dx = median(abs(diff(dop)));
    gfac = 2;
    if isfield(smooth, "lorentz_gamma_factor")
        gfac = smooth.lorentz_gamma_factor;
    end
    gamma = gfac * dx;
end

if isfinite(gamma) && gamma > 0
    chi_out = lorentz_smooth_1d_(dop, chi_in, gamma);
    info.enabled = true;
    info.method = "lorentz";
    info.gamma = gamma;
    fprintf("[chi smooth] method=lorentz gamma=%.6g doping units\n", gamma);
end

end

function plot_one_T_(R, par, plot_a1, plot_psi, plot_X)

fs = par.plot.fontSize;
xlo = min(R.dop_fwd);
xhi = max(R.dop_fwd);

if isfield(R, "B_T") && isfinite(R.B_T)
    title_info = sprintf('$B=%.3f\\,\\mathrm{T},\\; T=%.3f\\,\\mathrm{K},\\; q=(%d,%d)$', ...
        R.B_T, R.T, R.iq_pick, R.jq_pick);
else
    title_info = sprintf('$T=%.3f\\,\\mathrm{K},\\; q=(%d,%d)$', ...
        R.T, R.iq_pick, R.jq_pick);
end

line_f = '-';
line_b = '--';

if plot_a1
    figure('Color','w','Position',[80 80 720 560]);
    ax = axes; hold(ax,'on'); box(ax,'on');

    plot(ax, R.dop_fwd, R.a1_fwd, line_f, ...
        'Color','k', 'LineWidth',2.0);
    plot(ax, R.dop_bwd, R.a1_bwd, line_b, ...
        'Color','k', 'LineWidth',2.0);

    yline(ax, 0, '--', 'LineWidth',1.2, 'Color',[0.4 0.4 0.4]);

    xlim(ax,[xlo xhi]);
    grid(ax,'off');

    xlabel(ax,'doping $(10^{12}\,\mathrm{cm}^{-2})$','Interpreter','latex');
    ylabel(ax,'$a_1$','Interpreter','latex');
    title(ax,['$a_1$ vs doping, ' title_info], ...
        'Interpreter','latex','FontWeight','normal');

    set(ax, ...
        'FontSize',fs, ...
        'TickLabelInterpreter','latex', ...
        'LineWidth',1.3, ...
        'TickDir','in', ...
        'Box','on');
end

if plot_psi
    figure('Color','w','Position',[120 120 720 560]);
    ax = axes; hold(ax,'on'); box(ax,'on');

    plot(ax, R.dop_fwd, R.psi_f_abs_plot(:,1), line_f, ...
        'LineWidth',2.0);
    plot(ax, R.dop_bwd, R.psi_b_abs_plot(:,1), line_b, ...
        'LineWidth',2.0);

    xlim(ax,[xlo xhi]);
    grid(ax,'off');

    xlabel(ax,'doping $(10^{12}\,\mathrm{cm}^{-2})$','Interpreter','latex');
    ylabel(ax,'$|\psi|$','Interpreter','latex');
    title(ax,['$|\psi|$ hysteresis, ' title_info], ...
        'Interpreter','latex','FontWeight','normal');

    set(ax, ...
        'FontSize',fs, ...
        'TickLabelInterpreter','latex', ...
        'LineWidth',1.3, ...
        'TickDir','in', ...
        'Box','on');
end

if plot_X
    figure('Color','w','Position',[160 160 720 560]);
    ax = axes; hold(ax,'on'); box(ax,'on');

    plot(ax, R.dop_fwd, R.psi_f_abs_plot(:,2), line_f, ...
        'LineWidth',2.0);
    plot(ax, R.dop_bwd, R.psi_b_abs_plot(:,2), line_b, ...
        'LineWidth',2.0);

    xlim(ax,[xlo xhi]);
    grid(ax,'off');

    xlabel(ax,'doping $(10^{12}\,\mathrm{cm}^{-2})$','Interpreter','latex');
    ylabel(ax,'$|X|$','Interpreter','latex');
    title(ax,['$|X|$ hysteresis, ' title_info], ...
        'Interpreter','latex','FontWeight','normal');

    set(ax, ...
        'FontSize',fs, ...
        'TickLabelInterpreter','latex', ...
        'LineWidth',1.3, ...
        'TickDir','in', ...
        'Box','on');
end

end

function Tfolders = detect_T_folders_(root)

D = dir(root);
D = D([D.isdir]);

Ts = [];
Paths = strings(0,1);

for i = 1:numel(D)
    name = string(D(i).name);
    if name == "." || name == ".."
        continue;
    end

    tok = regexp(name, "^T([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)$", ...
        "tokens", "once");

    if isempty(tok)
        continue;
    end

    Tval = str2double(tok{1});
    if ~isfinite(Tval)
        continue;
    end

    Ts(end+1,1) = Tval; %#ok<AGROW>
    Paths(end+1,1) = string(fullfile(root, name)); %#ok<AGROW>
end

[Ts, ord] = sort(Ts);
Paths = Paths(ord);

Tfolders = struct('T', num2cell(Ts), 'path', cellstr(Paths));

end

function BTfolders = detect_BT_folders_(root, D_folder)

root = string(root);
D_folder = string(D_folder);

search_root = root;
if strlength(D_folder) > 0 && isfolder(fullfile(root, D_folder))
    search_root = string(fullfile(root, D_folder));
end

Bdirs = dir(fullfile(search_root, "B*T"));
Bdirs = Bdirs([Bdirs.isdir]);

Bvals = [];
Tvals = [];
Paths = strings(0,1);

for ib = 1:numel(Bdirs)
    bname = string(Bdirs(ib).name);
    Bval = parse_B_from_name_(bname);
    if ~isfinite(Bval)
        continue;
    end

    Tdirs = dir(fullfile(Bdirs(ib).folder, Bdirs(ib).name, "T*"));
    Tdirs = Tdirs([Tdirs.isdir]);

    for it = 1:numel(Tdirs)
        tname = string(Tdirs(it).name);
        Tval = parse_T_from_name_(tname);
        if ~isfinite(Tval)
            continue;
        end

        Bvals(end+1,1) = Bval; %#ok<AGROW>
        Tvals(end+1,1) = Tval; %#ok<AGROW>
        Paths(end+1,1) = string(fullfile(Tdirs(it).folder, Tdirs(it).name)); %#ok<AGROW>
    end
end

if isempty(Bvals)
    BTfolders = struct('B',{},'T',{},'path',{});
    return;
end

[~, ord] = sortrows([Bvals, Tvals], [1 2]);
Bvals = Bvals(ord);
Tvals = Tvals(ord);
Paths = Paths(ord);

BTfolders = struct('B', num2cell(Bvals), 'T', num2cell(Tvals), 'path', cellstr(Paths));

end

function Bval = parse_B_from_name_(name)

num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";
tok = regexp(string(name), "^B" + num + "T$", "tokens", "once", "ignorecase");
if isempty(tok)
    Bval = NaN;
else
    Bval = str2double(tok{1});
end

end

function Tval = parse_T_from_name_(name)

num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";
tok = regexp(string(name), "^T" + num + "$", "tokens", "once", "ignorecase");
if isempty(tok)
    Tval = NaN;
else
    Tval = str2double(tok{1});
end

end

function T_folder = find_T_folder_(Tfolders, T_target, T_tol)

Tall = [Tfolders.T].';

if isempty(Tall)
    T_folder = "";
    return;
end

[diffMin, id] = min(abs(Tall - T_target));

if isempty(id) || diffMin > max(T_tol, 1e-8)
    T_folder = "";
else
    T_folder = string(Tfolders(id).path);
end

end

function T_folder = find_BT_folder_(BTfolders, B_target, T_target, T_tol)

if isempty(BTfolders)
    T_folder = "";
    return;
end

Ball = [BTfolders.B].';
Tall = [BTfolders.T].';

B_tol = 1e-8;
ok = abs(Ball - B_target) <= B_tol;
if ~any(ok)
    T_folder = "";
    return;
end

idx = find(ok);
[diffMin, rel] = min(abs(Tall(idx) - T_target));

if isempty(rel) || diffMin > max(T_tol, 1e-8)
    T_folder = "";
else
    T_folder = string(BTfolders(idx(rel)).path);
end

end

function B = collect_mu_dop_chi_curve_(T_folder, iq_pick, jq_pick, chi_branch, broadening_meV)

if nargin < 4 || strlength(string(chi_branch)) == 0
    chi_branch = "average";
end
if nargin < 5
    broadening_meV = [];
end

file_pattern = chi_file_pattern_(chi_branch);
L = dir(fullfile(T_folder, "**", file_pattern));
L = L(~[L.isdir]);

if isempty(L)
    error("No %s under: %s", file_pattern, T_folder);
end

mu  = nan(numel(L),1);
dop = nan(numel(L),1);
chi = nan(numel(L),1);

n=0; dropH=0; dropTbl=0; dropQ=0; dropBroad=0; missingBroad=0;
useNameQ=0; useTableQ=0;

for k = 1:numel(L)
    fpath = string(fullfile(L(k).folder, L(k).name));
    folder = string(L(k).folder);

    [broad_file, has_broad_name] = parse_broadening_from_filename_(L(k).name);
    if ~isempty(broadening_meV)
        if ~has_broad_name
            missingBroad = missingBroad + 1;
            dropBroad = dropBroad + 1;
            continue;
        end
        if abs(broad_file - broadening_meV) > max(1e-9, 1e-6 * abs(broadening_meV))
            dropBroad = dropBroad + 1;
            continue;
        end
    end

    [iq_file, jq_file, has_q_name] = parse_q_from_filename_(L(k).name);
    if has_q_name && (iq_file ~= iq_pick || jq_file ~= jq_pick)
        dropQ = dropQ + 1;
        continue;
    end

    H = parse_header_mu_dop_(fpath);

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

    Cc = detect_cols_(M, fpath, iq_pick, jq_pick);

    if has_q_name
        id = find(isfinite(M(:, Cc.Re)), 1);
        if ~isempty(id)
            useNameQ = useNameQ + 1;
        end
    else
        id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
        if ~isempty(id)
            useTableQ = useTableQ + 1;
        end
    end

    if isempty(id)
        dropQ = dropQ + 1;
        continue;
    end

    n=n+1;
    mu(n)  = mu_use;
    dop(n) = H.dop;
    chi(n) = M(id, Cc.Re);
end

mu=mu(1:n);
dop=dop(1:n);
chi=chi(1:n);

ok = isfinite(mu) & isfinite(dop) & isfinite(chi);
mu=mu(ok);
dop=dop(ok);
chi=chi(ok);

fprintf("[collect] files=%d kept=%d | dropH=%d dropTbl=%d dropQ=%d dropBroad=%d missingBroad=%d | qName=%d qTable=%d\n", ...
    numel(L), numel(mu), dropH, dropTbl, dropQ, dropBroad, missingBroad, useNameQ, useTableQ);

if isempty(mu) && ~isempty(broadening_meV)
    available = list_broadening_values_(L);
    if isempty(available)
        error("No chi files with a broadening tag were found under %s.", T_folder);
    else
        error("No chi points left for broadening %.12g meV. Available broadenings in filenames: %s meV.", ...
            broadening_meV, strjoin(compose("%.12g", available), ", "));
    end
end

B = struct('mu',mu,'dop',dop,'chi',chi);

end

function values = list_broadening_values_(L)

values = [];
for ii = 1:numel(L)
    [b, ok] = parse_broadening_from_filename_(L(ii).name);
    if ok
        values(end+1) = b; %#ok<AGROW>
    end
end

values = unique(values(isfinite(values)));

end

function [broadening_meV, ok] = parse_broadening_from_filename_(fname)

broadening_meV = NaN;
ok = false;

tok = regexp(char(fname), ...
    'broadening([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)meV(?:_|\.|$)', ...
    'tokens', 'once');

if isempty(tok)
    return;
end

broadening_meV = str2double(tok{1});
ok = isfinite(broadening_meV);

end

function [iq, jq, ok] = parse_q_from_filename_(fname)

iq = NaN;
jq = NaN;
ok = false;

tok = regexp(char(fname), ...
    '_q([+\-]?\d+)_([+\-]?\d+)(?:_|\.|$)', ...
    'tokens', 'once');

if isempty(tok)
    return;
end

iq = str2double(tok{1});
jq = str2double(tok{2});
ok = isfinite(iq) && isfinite(jq);

if ok
    iq = round(iq);
    jq = round(jq);
end

end

function file_pattern = chi_file_pattern_(chi_branch)

chi_branch = string(chi_branch);

if chi_branch == "average"
    file_pattern = "chi_valley_*_spin_*.txt";
elseif startsWith(chi_branch, "chi_")
    file_pattern = chi_branch + "*.txt";
else
    file_pattern = "chi_" + chi_branch + "*.txt";
end

end

function H = parse_header_mu_dop_(fpath)

H = struct('mu',NaN,'dop',NaN);

fid = fopen(fpath,'r');
if fid < 0
    return;
end

c = onCleanup(@() fclose(fid)); %#ok<NASGU>
num  = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

for t = 1:3000
    ln = fgetl(fid);
    if ~ischar(ln)
        break;
    end

    s0 = strtrim(string(ln));
    if ~startsWith(s0,"#")
        break;
    end

    s = strtrim(erase(s0,"#"));

    if ~isfinite(H.dop)
        tok = regexp(s, "(?:^|\s)(?:doping|dop)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok)
            H.dop = str2double(tok{1});
        end
    end

    if ~isfinite(H.mu)
        tok = regexp(s, "(?:^|\s)(?:mu|EF|E_F)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok)
            H.mu = str2double(tok{1});
        end
    end
end

end

function mu = parse_mu_from_path_(folder)

folder = replace(string(folder),"\","/");
num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

tok = regexp(folder, "/mu\s*=?\s*" + num + "(?:/|$)", ...
    "tokens","once","ignorecase");

if ~isempty(tok)
    mu = str2double(tok{1});
else
    mu = NaN;
end

end

function C = detect_cols_(M, fpath, iq_pick, jq_pick)

C = detect_cols_from_header_(fpath);
if ~isempty(C)
    return;
end

candidates = [
    struct('iq',1,'jq',2,'Re',5)
    struct('iq',2,'jq',3,'Re',6)
];

for i = 1:numel(candidates)
    cand = candidates(i);
    if size(M,2) < max([cand.iq, cand.jq, cand.Re])
        continue;
    end

    if any(M(:,cand.iq)==iq_pick & M(:,cand.jq)==jq_pick)
        C = cand;
        return;
    end
end

if size(M,2) >= 8
    C = struct('iq',1,'jq',2,'Re',5);
elseif size(M,2) >= 6
    C = struct('iq',1,'jq',2,'Re',5);
else
    error("detect_cols_: unsupported chi table with %d columns: %s", ...
        size(M,2), fpath);
end

end

function C = detect_cols_from_header_(fpath)

C = [];

fid = fopen(fpath,'r');
if fid < 0
    return;
end

c = onCleanup(@() fclose(fid)); %#ok<NASGU>

for t = 1:3000
    ln = fgetl(fid);
    if ~ischar(ln)
        break;
    end

    s = strtrim(string(ln));
    if ~startsWith(s, "#")
        break;
    end

    s = strtrim(erase(s, "#"));

    if ~startsWith(lower(s), "iq ") && ...
       ~startsWith(lower(s), "idx ") && ...
       ~contains(lower(s), "chi_real") && ...
       ~contains(lower(s), "chi_re")
        continue;
    end

    toks = regexp(char(s), "\s+", "split");
    toks = toks(~cellfun(@isempty, toks));
    if isempty(toks)
        continue;
    end

    toks_norm = lower(string(toks));
    toks_norm = regexprep(toks_norm, "[^a-z0-9_]", "");
    if ~isempty(toks_norm) && toks_norm(1) == "columns"
        toks_norm = toks_norm(2:end);
    end

    iq_col = find(toks_norm == "iq", 1);
    jq_col = find(toks_norm == "jq", 1);

    re_names = ["chi_real", "chi_re", "chire", "real", "re"];
    re_col = [];
    for name = re_names
        re_col = find(toks_norm == name, 1);
        if ~isempty(re_col)
            break;
        end
    end

    if ~isempty(iq_col) && ~isempty(jq_col) && ~isempty(re_col)
        C = struct('iq',iq_col,'jq',jq_col,'Re',re_col);
        return;
    end
end

end

function [mu2, chi2, dop2] = collapse_same_mu_(mu, chi, dop)

mu  = double(mu(:));
chi = double(chi(:));
dop = double(dop(:));

ok = isfinite(mu) & isfinite(chi) & isfinite(dop);
mu=mu(ok);
chi=chi(ok);
dop=dop(ok);

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
dop=dop(ok);
chi=chi(ok);
mu=mu(ok);

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

if ~isfinite(scan_start)
    scan_start = dmin_all;
end
if ~isfinite(scan_end)
    scan_end = dmax_all;
end

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
    N = double(hyst.N);
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

function y_s = lorentz_smooth_1d_(x, y, gamma)

x = x(:);
y = y(:);
N = numel(x);
y_s = zeros(N,1);

for i = 1:N
    dx = x(i) - x;
    w  = (gamma/pi) ./ (dx.^2 + gamma^2);
    sw = sum(w);

    if sw > 0
        y_s(i) = sum(w .* y) / sw;
    else
        y_s(i) = y(i);
    end
end

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

if R < 3 || C < 3
    return;
end

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

function write_hyst_data_BT(R, ib, it, nB, nT, output_dir, run_stamp)

% ============================================================
% Filename rule:
% nB = 1, nT = 1  -> hyst_DBT_<timestamp>.mat
% nB > 1          -> append _allB
% nT > 1          -> append _allT
% ============================================================

suffix = "";

if nB > 1
    suffix = suffix + "_allB";
end

if nT > 1
    suffix = suffix + "_allT";
end

fname = fullfile(output_dir, "hyst_DBT" + suffix + "_" + run_stamp + ".mat");

entry = struct();

entry.T = R.T;
entry.T_folder = R.T_folder;
entry.B_T = R.B_T;
entry.B_folder = R.B_folder;
entry.D_folder = R.D_folder;
entry.chi_branch = R.chi_branch;
entry.iq_pick = R.iq_pick;
entry.jq_pick = R.jq_pick;
entry.dmu_eV = R.dmu_eV;
entry.doping_scale = R.doping_scale;
entry.broadening_meV = R.broadening_meV;
entry.par = R.par;
entry.run_params = make_run_params_record_(R, R.par);

entry.dop_fwd = R.dop_fwd;
entry.dop_bwd = R.dop_bwd;

entry.psi_fwd = R.psi_fwd(:,1);
entry.psi_bwd = R.psi_bwd(:,1);

entry.X_fwd = R.psi_fwd(:,2);
entry.X_bwd = R.psi_bwd(:,2);

entry.abspsi_fwd = abs(R.psi_fwd(:,1));
entry.abspsi_bwd = abs(R.psi_bwd(:,1));

entry.absX_fwd = abs(R.psi_fwd(:,2));
entry.absX_bwd = abs(R.psi_bwd(:,2));

entry.psi_f_plot = R.psi_f_abs_plot(:,1);
entry.psi_b_plot = R.psi_b_abs_plot(:,1);

entry.X_f_plot = R.psi_f_abs_plot(:,2);
entry.X_b_plot = R.psi_b_abs_plot(:,2);

entry.a1_fwd = R.a1_fwd;
entry.a1_bwd = R.a1_bwd;

entry.chi_fwd = R.chi_fwd;
entry.chi_bwd = R.chi_bwd;

entry.mu_raw = R.mu_raw;
entry.dop_raw_from_file = R.dop_raw_from_file;
entry.dop_raw = R.dop_raw;
entry.chi_raw = R.chi_raw;

entry.mu_used = R.mu_used;
entry.dop_used = R.dop_used;
entry.chi_used = R.chi_used;
entry.chi_unsmoothed = R.chi_unsmoothed;
entry.chi_smooth_info = R.chi_smooth_info;
entry.gamma_used = R.gamma_used;

if isfile(fname)
    S = load(fname, "HYST");
    HYST = S.HYST;
    sz = size(HYST);
    if numel(sz) < 2
        sz(end+1:2) = 1;
    end
    field_mismatch = ~isequal(sort(fieldnames(HYST)), sort(fieldnames(entry)));
    size_mismatch = any(sz(1:2) < [nB, nT]);
    if field_mismatch || size_mismatch
        HYST = repmat(entry, nB, nT);
    end
else
    HYST = repmat(entry, nB, nT);
end

HYST(ib, it) = entry;

RUN = entry.run_params;
save(fname, "HYST", "RUN");

fprintf("[write MAT] %s | B=%.3f T | T=%.3f K\n", ...
    fname, R.B_T, R.T);

end
