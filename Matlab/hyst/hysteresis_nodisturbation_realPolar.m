function H = hysteresis_nodisturbation_realPolar()
clc; close all;

% =========================
% USER SETTINGS
% =========================
T_target = 6.5;              % no T folder, but still needed for coef.eval
polar_list_meV = [];         % [] = use all polar folders

plot_a1  = true;
plot_psi = true;
plot_X   = true;

default_root = "E:/rg_master/data/";

% =========================
% Init parameters
% =========================
par  = make_realchi_params(true);
coef = make_realchi_coeff(par);

iq_pick = round(par.iq_pick);
jq_pick = round(par.jq_pick);

% =========================
% Select root folder
% root / polarXXX / muXXX / chi*.txt
% =========================
if ~isfolder(default_root)
    default_root = string(pwd);
end

root_dir = uigetdir(default_root, ...
    "Select root folder containing polarXXX/muXXX/chi*.txt");

if isequal(root_dir,0)
    H = [];
    return;
end

root_dir = string(root_dir);

polar_folders = find_polar_folders(root_dir);

if isempty(polar_folders)
    error("No polar folders found under: %s", root_dir);
end

all_polar = [polar_folders.polar_meV];

if isempty(polar_list_meV)
    polar_use = all_polar;
else
    polar_use = polar_list_meV(:).';
end

fprintf("[root] %s\n", root_dir);
fprintf("[q] (%d,%d)\n", iq_pick, jq_pick);
fprintf("[T] %.6g K\n", T_target);
fprintf("[polar used] ");
fprintf("%.6g ", polar_use);
fprintf("\n");

results = struct([]);

for ip = 1:numel(polar_use)

    polar_meV = polar_use(ip);
    polar_dir = match_polar_folder(polar_folders, polar_meV);

    if strlength(polar_dir) == 0
        warning("Skip polar %.6g meV: folder not found.", polar_meV);
        continue;
    end

    fprintf("\n================ polar = %.6g meV ================\n", polar_meV);
    fprintf("folder: %s\n", polar_dir);

    try
        R = run_hysteresis_for_one_polar( ...
            polar_dir, polar_meV, T_target, iq_pick, jq_pick, par, coef);
    catch ME
        warning("Skip polar %.6g meV: %s", polar_meV, ME.message);
        continue;
    end

    results(end+1) = R; %#ok<AGROW>

    save_hysteresis_result(R);
    plot_hysteresis_result(R, par, plot_a1, plot_psi, plot_X);
end

H = struct();
H.root_dir = root_dir;
H.T_target = T_target;
H.iq_pick = iq_pick;
H.jq_pick = jq_pick;
H.polar_list_meV = polar_use;
H.results = results;
H.par = par;

end

% ============================================================
% Run hysteresis for one real polar folder
% ============================================================
function R = run_hysteresis_for_one_polar(polar_dir, polar_meV, T_target, iq_pick, jq_pick, par, coef)

D = collect_chi_curve_from_polar_folder(polar_dir, iq_pick, jq_pick);

mu_raw  = D.mu(:);
dop_raw = D.doping(:);
chi_raw = D.chi(:);

if numel(mu_raw) < 5
    error("Too few valid chi points.");
end

fprintf("[raw] points=%d | mu=[%.6g, %.6g] | doping=[%.6g, %.6g]\n", ...
    numel(mu_raw), min(mu_raw), max(mu_raw), min(dop_raw), max(dop_raw));

[dop_used, ord] = sort(dop_raw);
chi_used = chi_raw(ord);
mu_used  = mu_raw(ord);

[dop_used, chi_used, mu_used] = average_duplicate_doping(dop_used, chi_used, mu_used);

chi_of_dop = @(d) interp1(dop_used, chi_used, d, "linear", "extrap");

[dop_fwd, dop_bwd] = build_hyst_doping_grid_(dop_used, par.hyst);
N = numel(dop_fwd);

fprintf("[scan] fwd=[%.6g, %.6g], bwd=[%.6g, %.6g], N=%d\n", ...
    dop_fwd(1), dop_fwd(end), dop_bwd(1), dop_bwd(end), N);

scan = make_scan_options(par);
opt  = make_min_options(par);

psi_f = nan(N,2);
psi_b = nan(N,2);

a1_f  = nan(N,1);
a1_b  = nan(N,1);

chi_f = nan(N,1);
chi_b = nan(N,1);

% =========================
% Forward scan
% =========================
prev = [par.hyst.psi1_0, par.hyst.psi2_0];

for i = 1:N

    dop = dop_fwd(i);
    chi_now = chi_of_dop(dop);
    C = coef.eval(T_target, dop, chi_now);

    pick = choose_local_minimum(C, par, scan, opt, prev, true);

    psi_f(i,:) = pick;
    prev = pick;

    a1_f(i) = C.a1;
    chi_f(i) = chi_now;
end

% =========================
% Backward scan
% =========================
prev = psi_f(end,:);

for i = 1:N

    dop = dop_bwd(i);
    chi_now = chi_of_dop(dop);
    C = coef.eval(T_target, dop, chi_now);

    pick = choose_local_minimum(C, par, scan, opt, prev, false);

    psi_b(i,:) = pick;
    prev = pick;

    a1_b(i) = C.a1;
    chi_b(i) = chi_now;
end

psi_f_abs = abs(psi_f);
psi_b_abs = abs(psi_b);

[psi_f_plot, psi_b_plot, gamma_used] = apply_lorentz_on_doping_( ...
    dop_fwd, psi_f_abs, dop_bwd, psi_b_abs, par.smooth);

R = struct();

R.T = T_target;
R.polar = polar_meV;
R.polar_meV = polar_meV;
R.real_polar = polar_meV;

R.iq_pick = iq_pick;
R.jq_pick = jq_pick;
R.polar_dir = polar_dir;

R.mu_raw = mu_raw;
R.dop_raw = dop_raw;
R.chi_raw = chi_raw;

R.mu_used = mu_used;
R.dop_used = dop_used;
R.chi_used = chi_used;

R.dop_fwd = dop_fwd;
R.dop_bwd = dop_bwd;

R.psi_fwd = psi_f;
R.psi_bwd = psi_b;

R.psi_f_plot = psi_f_plot(:,1);
R.psi_b_plot = psi_b_plot(:,1);

R.X_f_plot = psi_f_plot(:,2);
R.X_b_plot = psi_b_plot(:,2);

R.a1_fwd = a1_f;
R.a1_bwd = a1_b;

R.chi_fwd = chi_f;
R.chi_bwd = chi_b;

R.gamma_used = gamma_used;

end

% ============================================================
% Collect chi from:
% polar_dir / muXXX / chi*.txt
% ============================================================
function D = collect_chi_curve_from_polar_folder(polar_dir, iq_pick, jq_pick)

files = dir(fullfile(polar_dir, "mu*", "chi*.txt"));

if isempty(files)
    error("No chi*.txt found under %s/mu*/", polar_dir);
end

mu = [];
doping = [];
chi = [];

drop_header = 0;
drop_table  = 0;
drop_q      = 0;

for k = 1:numel(files)

    fpath = string(fullfile(files(k).folder, files(k).name));

    meta = parse_header_mu_doping(fpath);
    mu_val = parse_mu_from_folder(files(k).folder);

    if ~isfinite(mu_val)
        mu_val = meta.mu;
    end

    if ~isfinite(mu_val) || ~isfinite(meta.doping)
        drop_header = drop_header + 1;
        continue;
    end

    M = read_chi_table(fpath);

    if isempty(M) || size(M,2) < 6
        drop_table = drop_table + 1;
        continue;
    end

    C = detect_chi_columns(M);
    id = find(M(:,C.iq) == iq_pick & M(:,C.jq) == jq_pick, 1);

    if isempty(id)
        drop_q = drop_q + 1;
        continue;
    end

    mu(end+1,1) = mu_val; %#ok<AGROW>
    doping(end+1,1) = meta.doping; %#ok<AGROW>
    chi(end+1,1) = M(id,C.Re); %#ok<AGROW>
end

ok = isfinite(mu) & isfinite(doping) & isfinite(chi);

D = struct();
D.mu = mu(ok);
D.doping = doping(ok);
D.chi = chi(ok);

fprintf("[collect] files=%d kept=%d | dropHeader=%d dropTable=%d dropQ=%d\n", ...
    numel(files), numel(D.mu), drop_header, drop_table, drop_q);

end

% ============================================================
% Find polar folders
% Accept:
% polar0
% polar0.5
% polar_meV0.5
% polar_0.5
% ============================================================
function polar_folders = find_polar_folders(root_dir)

items = dir(root_dir);
items = items([items.isdir]);

polar_vals = [];
paths = strings(0,1);

num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

for i = 1:numel(items)

    name = string(items(i).name);

    if name == "." || name == ".."
        continue;
    end

    tok = regexp(name, "^polar(?:_meV)?_?" + num + "$", "tokens", "once");

    if isempty(tok)
        continue;
    end

    p = str2double(tok{1});

    if ~isfinite(p)
        continue;
    end

    polar_vals(end+1,1) = p; %#ok<AGROW>
    paths(end+1,1) = string(fullfile(items(i).folder, items(i).name)); %#ok<AGROW>
end

[polar_vals, ord] = sort(polar_vals);
paths = paths(ord);

polar_folders = struct( ...
    "polar_meV", num2cell(polar_vals), ...
    "path", cellstr(paths));

end

function polar_dir = match_polar_folder(polar_folders, polar_target)

all_polar = [polar_folders.polar_meV];

[diff_min, id] = min(abs(all_polar - polar_target));

if isempty(id) || diff_min > 1e-10
    polar_dir = "";
else
    polar_dir = string(polar_folders(id).path);
end

end

% ============================================================
% Choose local minimum
% ============================================================
function pick = choose_local_minimum(C, par, scan, opt, prev, allow_zero_if_a1_positive)

if allow_zero_if_a1_positive && C.a1 > 0
    pick = [0,0];
    return;
end

F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);

if isempty(mins)
    [p1,p2,~] = minimize_free_energy( ...
        C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda, ...
        prev(1), prev(2), false);
    pick = [p1,p2];
else
    pick = pick_closest_min_(mins, prev);
end

end

function scan = make_scan_options(par)

scan = struct();
scan.Npsi1 = par.min.Npsi1;
scan.Npsi2 = par.min.Npsi2;
scan.use_gaussian_smooth = par.min.use_gaussian_smooth;
scan.smooth_sigma = par.min.smooth_sigma;
scan.keep_topK = par.min.keep_topK;

end

function opt = make_min_options(par)

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

% ============================================================
% Save MAT
% ============================================================
function save_hysteresis_result(R)

fname = fullfile(pwd, "hyst_allRealPolar.mat");

entry = struct();

entry.T = R.T;
entry.polar = R.polar;
entry.polar_meV = R.polar_meV;
entry.real_polar = R.real_polar;

entry.artificial_polar = NaN;
entry.dmu_eV = NaN;

entry.iq_pick = R.iq_pick;
entry.jq_pick = R.jq_pick;
entry.polar_dir = R.polar_dir;

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

entry.psi_f_plot = R.psi_f_plot;
entry.psi_b_plot = R.psi_b_plot;

entry.X_f_plot = R.X_f_plot;
entry.X_b_plot = R.X_b_plot;

entry.a1_fwd = R.a1_fwd;
entry.a1_bwd = R.a1_bwd;

entry.chi_fwd = R.chi_fwd;
entry.chi_bwd = R.chi_bwd;

entry.mu_raw = R.mu_raw;
entry.dop_raw = R.dop_raw;
entry.chi_raw = R.chi_raw;

entry.mu_used = R.mu_used;
entry.dop_used = R.dop_used;
entry.chi_used = R.chi_used;

if isfile(fname)
    S = load(fname, "HYST");
    HYST = S.HYST;
    HYST(end+1) = entry;
else
    HYST = entry;
end

save(fname, "HYST");

fprintf("[write] %s | polar=%.6g meV\n", fname, R.polar);

end

% ============================================================
% Plot
% ============================================================
function plot_hysteresis_result(R, par, plot_a1, plot_psi, plot_X)

fs = par.plot.fontSize;

xlo = min([R.dop_fwd(:); R.dop_bwd(:)]);
xhi = max([R.dop_fwd(:); R.dop_bwd(:)]);

title_info = sprintf('$T=%.3f\\,\\mathrm{K},\\; polar=%.3f\\,\\mathrm{meV},\\; q=(%d,%d)$', ...
    R.T, R.polar, R.iq_pick, R.jq_pick);

if plot_a1
    figure('Color','w','Position',[80 80 720 560]);
    ax = axes; hold(ax,'on'); box(ax,'on');

    plot(ax, R.dop_fwd, R.a1_fwd, '-', 'Color','k', 'LineWidth',2.0);
    plot(ax, R.dop_bwd, R.a1_bwd, '--', 'Color','k', 'LineWidth',2.0);
    yline(ax, 0, '--', 'LineWidth',1.2, 'Color',[0.4 0.4 0.4]);

    xlim(ax,[xlo xhi]);
    xlabel(ax,'doping $(10^{12}\,\mathrm{cm}^{-2})$','Interpreter','latex');
    ylabel(ax,'$a_1$','Interpreter','latex');
    title(ax,['$a_1$ vs doping, ' title_info], ...
        'Interpreter','latex','FontWeight','normal');

    set(ax,'FontSize',fs,'TickLabelInterpreter','latex', ...
        'LineWidth',1.3,'TickDir','in','Box','on');
end

if plot_psi
    figure('Color','w','Position',[120 120 720 560]);
    ax = axes; hold(ax,'on'); box(ax,'on');

    plot(ax, R.dop_fwd, R.psi_f_plot, '-', 'LineWidth',2.0);
    plot(ax, R.dop_bwd, R.psi_b_plot, '--', 'LineWidth',2.0);

    xlim(ax,[xlo xhi]);
    xlabel(ax,'doping $(10^{12}\,\mathrm{cm}^{-2})$','Interpreter','latex');
    ylabel(ax,'$|\psi|$','Interpreter','latex');
    title(ax,['$|\psi|$ hysteresis, ' title_info], ...
        'Interpreter','latex','FontWeight','normal');

    set(ax,'FontSize',fs,'TickLabelInterpreter','latex', ...
        'LineWidth',1.3,'TickDir','in','Box','on');
end

if plot_X
    figure('Color','w','Position',[160 160 720 560]);
    ax = axes; hold(ax,'on'); box(ax,'on');

    plot(ax, R.dop_fwd, R.X_f_plot, '-', 'LineWidth',2.0);
    plot(ax, R.dop_bwd, R.X_b_plot, '--', 'LineWidth',2.0);

    xlim(ax,[xlo xhi]);
    xlabel(ax,'doping $(10^{12}\,\mathrm{cm}^{-2})$','Interpreter','latex');
    ylabel(ax,'$|X|$','Interpreter','latex');
    title(ax,['$|X|$ hysteresis, ' title_info], ...
        'Interpreter','latex','FontWeight','normal');

    set(ax,'FontSize',fs,'TickLabelInterpreter','latex', ...
        'LineWidth',1.3,'TickDir','in','Box','on');
end

end

% ============================================================
% Parsers
% ============================================================
function meta = parse_header_mu_doping(fpath)

meta = struct("mu",NaN,"doping",NaN);

fid = fopen(fpath,"r");
if fid < 0
    return;
end

c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

for i = 1:4000

    ln = fgetl(fid);

    if ~ischar(ln)
        break;
    end

    s0 = strtrim(string(ln));

    if strlength(s0) == 0
        continue;
    end

    if ~startsWith(s0,"#")
        break;
    end

    s = strtrim(erase(s0,"#"));

    if ~isfinite(meta.doping)
        tok = regexp(s, "(?:^|\s)(?:doping|dop)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok)
            meta.doping = str2double(tok{1});
        end
    end

    if ~isfinite(meta.mu)
        tok = regexp(s, "(?:^|\s)(?:mu|EF|E_F)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok)
            meta.mu = str2double(tok{1});
        end
    end
end

end

function mu = parse_mu_from_folder(folder)

folder = replace(string(folder), "\", "/");

num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

tok = regexp(folder, "/mu_?" + num + "(?:/|$)", "tokens","once","ignorecase");

if isempty(tok)
    tok = regexp(folder, "/mu=?" + num + "(?:/|$)", "tokens","once","ignorecase");
end

if isempty(tok)
    mu = NaN;
else
    mu = str2double(tok{1});
end

end

function M = read_chi_table(fpath)

try
    M = readmatrix(fpath, "FileType","text", "CommentStyle","#");
catch
    M = [];
end

end

function C = detect_chi_columns(M)

if size(M,2) >= 8
    C = struct("iq",2,"jq",3,"Re",6);
else
    C = struct("iq",1,"jq",2,"Re",5);
end

end

% ============================================================
% Average duplicate doping
% ============================================================
function [dop2, chi2, mu2] = average_duplicate_doping(dop, chi, mu)

dop = double(dop(:));
chi = double(chi(:));
mu  = double(mu(:));

ok = isfinite(dop) & isfinite(chi) & isfinite(mu);

dop = dop(ok);
chi = chi(ok);
mu  = mu(ok);

dop_key = round(dop, 12);

[dop_u,~,ic] = unique(dop_key);

chi_m = accumarray(ic, chi, [], @mean);
mu_m  = accumarray(ic, mu,  [], @mean);

[dop2,ord] = sort(dop_u);
chi2 = chi_m(ord);
mu2  = mu_m(ord);

end

% ============================================================
% Build scan grid
% ============================================================
function [dop_fwd, dop_bwd] = build_hyst_doping_grid_(dop0, hyst)

dop0 = double(dop0(:));
dop0 = dop0(isfinite(dop0));

if isempty(dop0)
    error("dop0 is empty.");
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

dlo = max(dmin_all, min(scan_start, scan_end));
dhi = min(dmax_all, max(scan_start, scan_end));

if ~(isfinite(dlo) && isfinite(dhi) && dhi > dlo)
    error("Requested scan window has no overlap with data range.");
end

if mode == "step"

    step = double(hyst.step);

    if ~isfinite(step) || step <= 0
        error("hyst.step must be positive.");
    end

    dop_base = (dlo:step:dhi).';

    if abs(dop_base(end) - dhi) > 1e-12
        dop_base(end+1,1) = dhi;
    end

elseif mode == "N"

    N = double(hyst.N);

    if ~isfinite(N) || N < 2
        error("hyst.N must be >= 2.");
    end

    dop_base = linspace(dlo, dhi, round(N)).';

else
    error("Unsupported hyst.grid_mode.");
end

if direction == "ascend"
    dop_fwd = dop_base;
else
    dop_fwd = flipud(dop_base);
end

dop_bwd = flipud(dop_fwd);

end

% ============================================================
% Smoothing
% ============================================================
function [Yf_plot, Yb_plot, gamma_used] = apply_lorentz_on_doping_(xf, Yf, xb, Yb, smooth)

Yf_plot = Yf;
Yb_plot = Yb;
gamma_used = NaN;

if ~isfield(smooth,'use_lorentz') || ~smooth.use_lorentz
    return;
end

dx = median(abs(diff(xf)));

if ~isfinite(dx) || dx <= 0
    return;
end

gfac = 0;

if isfield(smooth,'lorentz_gamma_factor')
    gfac = smooth.lorentz_gamma_factor;
end

if ~isfinite(gfac) || gfac < 0
    return;
end

if gfac == 0
    gamma = 1.5 * dx;
else
    gamma = gfac * dx;
end

gamma_used = gamma;

Yf_plot = zeros(size(Yf));
Yb_plot = zeros(size(Yb));

for c = 1:size(Yf,2)
    Yf_plot(:,c) = lorentz_smooth_1d_(xf, Yf(:,c), gamma);
end

for c = 1:size(Yb,2)
    Yb_plot(:,c) = lorentz_smooth_1d_(xb, Yb(:,c), gamma);
end

end

function y_s = lorentz_smooth_1d_(x, y, gamma)

x = x(:);
y = y(:);

N = numel(x);
y_s = zeros(N,1);

for i = 1:N
    dx = x(i) - x;
    w = (gamma/pi) ./ (dx.^2 + gamma^2);
    y_s(i) = sum(w .* y) / sum(w);
end

end

% ============================================================
% Local minima utilities
% ============================================================
function pick = pick_closest_min_(mins, prev)

d = sqrt(sum((mins - prev).^2, 2));
[~,id] = min(d);
pick = mins(id,:);

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
    seeds = [seeds; 0,0; d,0; -d,0; 0,d; 0,-d];
end

obj = @(x) safe_F_(F, x(1), x(2), psi1_lim, psi2_lim);

fopt = optimset('Display','off', ...
    'MaxIter', opt.max_iter, ...
    'TolFun', opt.tol_fun, ...
    'TolX', opt.tol_x);

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
    Hh = hessian_fd_(@(u,v) F(u,v), x(1), x(2), opt.fd_h);
    ev = eig((Hh + Hh')/2);

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

fxx = (f(x+h,y) - 2*f00 + f(x-h,y)) / h^2;
fyy = (f(x,y+h) - 2*f00 + f(x,y-h)) / h^2;

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
n2 = A(3:R,   2:C-1);
n3 = A(2:R-1, 1:C-2);
n4 = A(2:R-1, 3:C);

n5 = A(1:R-2, 1:C-2);
n6 = A(1:R-2, 3:C);
n7 = A(3:R,   1:C-2);
n8 = A(3:R,   3:C);

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
x = -rad:rad;

g = exp(-(x.^2)/(2*sigma^2));
g = g / sum(g);

B = conv2(A, g, 'same');
B = conv2(B, g', 'same');

end