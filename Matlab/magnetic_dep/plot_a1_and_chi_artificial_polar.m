function H = plot_a1_and_chi_artificial_polar()

clc; close all;

% ============================================================
% USER SETTINGS
% ============================================================
D_root = "E:\rg_master\data";

T_target = 6.5;
T_tol    = 1e-6;

iq_pick = 3;
jq_pick = -3;

x_mode = "doping";          % "mu" or "doping"
use_abs_chi = false;

mu_range  = [];
dop_range = [-2, 0];

polar_list_meV = 0:0.5:6;  % artificial polar list

mu_unit_in_header = "eV";  % "eV" or "meV"

FS = 15;
LW = 2.0;
use_colorbar = true;

% ============================================================
% Coefficients for a1
% ============================================================
par  = make_realchi_params(true);
coef = make_realchi_coeff(par);

% ============================================================
% Select baseline folder
% ============================================================
if ~isfolder(D_root)
    D_root = pwd;
end

baseline_dir = uigetdir(char(D_root), ...
    "Select BASELINE polar=0 folder containing chi*.txt");

if isequal(baseline_dir,0)
    return;
end

baseline_dir = string(baseline_dir);

% ============================================================
% Collect baseline curve
% ============================================================
[mu0, dop0, chi0] = collect_curve_from_folder( ...
    baseline_dir, T_target, T_tol, iq_pick, jq_pick, mu_range, dop_range);

if numel(mu0) < 3
    error("Baseline curve has too few points.");
end

mu0  = double(mu0(:));
dop0 = double(dop0(:));
chi0 = double(chi0(:));

[mu0, ord] = sort(mu0);
dop0 = dop0(ord);
chi0 = chi0(ord);

Fchi = griddedInterpolant(mu0, chi0, "linear", "none");

switch lower(string(mu_unit_in_header))
    case "ev"
        polar_to_mu = 1e-3;
    case "mev"
        polar_to_mu = 1.0;
    otherwise
        error('mu_unit_in_header must be "eV" or "meV".');
end

polar_use = unique(double(polar_list_meV(:)));
polar_use = polar_use(isfinite(polar_use));
polar_use = sort(polar_use);

pmin = min(polar_use);
pmax = max(polar_use);
if pmin == pmax
    pmax = pmin + 1e-12;
end

cmap = blue_gradient(256);
get_color = @(p) polar_to_color_(p, pmin, pmax, cmap);

% ============================================================
% Prepare output
% ============================================================
H = struct();
H.baseline_dir = baseline_dir;
H.T_target = T_target;
H.T_tol = T_tol;
H.iq_pick = iq_pick;
H.jq_pick = jq_pick;
H.x_mode = x_mode;
H.polar_list_meV = polar_use;
H.curves = struct([]);

% ============================================================
% Figure 1: chi
% ============================================================
fig_chi = figure("Color","w","Units","pixels","Position",[120 120 900 560]);
ax_chi = axes(fig_chi);
hold(ax_chi,"on");
box(ax_chi,"on");
grid(ax_chi,"off");

set(ax_chi, ...
    "FontSize",FS, ...
    "LineWidth",1.3, ...
    "TickDir","in", ...
    "Box","on");

% ============================================================
% Figure 2: a1
% ============================================================
fig_a1 = figure("Color","w","Units","pixels","Position",[160 140 900 560]);
ax_a1 = axes(fig_a1);
hold(ax_a1,"on");
box(ax_a1,"on");
grid(ax_a1,"off");

set(ax_a1, ...
    "FontSize",FS, ...
    "LineWidth",1.3, ...
    "TickDir","in", ...
    "Box","on");

% ============================================================
% Loop over artificial polar
% ============================================================
for ip = 1:numel(polar_use)

    p = polar_use(ip);
    dmu = 0.5 * p * polar_to_mu;

    chi_plus  = Fchi(mu0 + dmu);
    chi_minus = Fchi(mu0 - dmu);
    chi_art   = 0.5 * (chi_plus + chi_minus);

    ok = isfinite(mu0) & isfinite(dop0) & isfinite(chi_art);

    if nnz(ok) < 3
        fprintf("[skip] artificial polar %.6g meV has too few valid points\n", p);
        continue;
    end

    mu_plot  = mu0(ok);
    dop_plot = dop0(ok);
    chi_plot = chi_art(ok);

    if use_abs_chi
        chi_for_plot = abs(chi_plot);
    else
        chi_for_plot = chi_plot;
    end

    if x_mode == "mu"
        x_plot = mu_plot;
    else
        x_plot = dop_plot;
    end

    % average duplicate x
    [x_plot, chi_for_plot, dop_plot, mu_plot] = average_by_x4_( ...
        x_plot, chi_for_plot, dop_plot, mu_plot);

    % compute a1 from artificial chi
    a1_plot = nan(size(x_plot));
    for i = 1:numel(x_plot)
        C = coef.eval(T_target, dop_plot(i), chi_for_plot(i));
        a1_plot(i) = C.a1;
    end

    [x_plot, ord] = sort(x_plot);
    chi_for_plot = chi_for_plot(ord);
    a1_plot = a1_plot(ord);
    dop_plot = dop_plot(ord);
    mu_plot = mu_plot(ord);

    if abs(p) < 1e-12
        line_style = "--";
        this_color = [0.55 0.55 0.55];
        this_lw = 2.4;
    else
        line_style = "-";
        this_color = get_color(p);
        this_lw = LW;
    end

    plot(ax_chi, x_plot, chi_for_plot, line_style, ...
        "LineWidth", this_lw, ...
        "Color", this_color);

    plot(ax_a1, x_plot, a1_plot, line_style, ...
        "LineWidth", this_lw, ...
        "Color", this_color);

    H.curves(end+1).polar_meV = p; %#ok<AGROW>
    H.curves(end).x = x_plot;
    H.curves(end).chi = chi_for_plot;
    H.curves(end).a1 = a1_plot;
    H.curves(end).doping = dop_plot;
    H.curves(end).mu = mu_plot;
end

% ============================================================
% Labels
% ============================================================
if x_mode == "mu"
    xlab = "\mu (eV)";
else
    xlab = "doping (10^{12} cm^{-2})";
end

xlabel(ax_chi, xlab, "Interpreter","tex", "FontSize",FS);
xlabel(ax_a1,  xlab, "Interpreter","tex", "FontSize",FS);

if use_abs_chi
    ylabel(ax_chi, "|Re(\chi)|", "Interpreter","tex", "FontSize",FS);
else
    ylabel(ax_chi, "Re(\chi)", "Interpreter","tex", "FontSize",FS);
end

ylabel(ax_a1, "$a_1$", "Interpreter","latex", "FontSize",FS);

title(ax_chi, sprintf("Artificial $\\chi$, $T=%.4g$ K, $q=(%d,%d)$", ...
    T_target, iq_pick, jq_pick), ...
    "Interpreter","latex", "FontWeight","normal", "FontSize",FS);

title(ax_a1, sprintf("$a_1$ from artificial $\\chi$, $T=%.4g$ K, $q=(%d,%d)$", ...
    T_target, iq_pick, jq_pick), ...
    "Interpreter","latex", "FontWeight","normal", "FontSize",FS);

yline(ax_a1, 0, "--", "LineWidth",1.1);

% ============================================================
% Colorbars
% ============================================================
if use_colorbar
    colormap(ax_chi, cmap);
    caxis(ax_chi, [pmin pmax]);
    cb1 = colorbar(ax_chi);
    cb1.Label.String = "artificial polar (meV)";
    cb1.TickDirection = "in";
    cb1.LineWidth = 1.2;

    colormap(ax_a1, cmap);
    caxis(ax_a1, [pmin pmax]);
    cb2 = colorbar(ax_a1);
    cb2.Label.String = "artificial polar (meV)";
    cb2.TickDirection = "in";
    cb2.LineWidth = 1.2;
end

H.fig_chi = fig_chi;
H.ax_chi = ax_chi;
H.fig_a1 = fig_a1;
H.ax_a1 = ax_a1;

end

% ============================================================
% collect baseline chi curve
% ============================================================
function [mu0, dop0, chi0] = collect_curve_from_folder(folder, T_target, T_tol, ...
    iq_pick, jq_pick, mu_range, dop_range)

files0 = dir(fullfile(folder, "**", "chi*.txt"));
files0 = files0(~[files0.isdir]);

if isempty(files0)
    error("No chi*.txt found under selected folder: %s", folder);
end

T_lo = T_target - T_tol;
T_hi = T_target + T_tol;

mu_list  = [];
dop_list = [];
chi_list = [];

for k = 1:numel(files0)

    fpath = string(fullfile(files0(k).folder, files0(k).name));
    meta = parse_header_T_doping_mu_(fpath);

    if ~isfinite(meta.T), continue; end
    if ~(meta.T >= T_lo && meta.T <= T_hi), continue; end
    if ~pass_range_AND_(meta.mu, mu_range), continue; end
    if ~pass_range_AND_(meta.doping, dop_range), continue; end

    M = read_numeric_skiphash_(fpath);
    if isempty(M) || size(M,2) < 6
        continue;
    end

    C = detect_cols_chi_(M);
    id = find(M(:,C.iq) == iq_pick & M(:,C.jq) == jq_pick, 1);

    if isempty(id)
        continue;
    end

    mu_list(end+1,1)  = meta.mu; %#ok<AGROW>
    dop_list(end+1,1) = meta.doping; %#ok<AGROW>
    chi_list(end+1,1) = M(id, C.Re); %#ok<AGROW>
end

[mu0, dop0, chi0] = collapse_same_mu_(mu_list, dop_list, chi_list);

end

% ============================================================
% parsers
% ============================================================
function meta = parse_header_T_doping_mu_(fpath)

meta = struct("T",NaN, "doping",NaN, "mu",NaN);

fid = fopen(fpath,"r");
if fid < 0
    return;
end

c = onCleanup(@() fclose(fid)); %#ok<NASGU>
num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';

for t = 1:4000

    ln = fgetl(fid);
    if ~ischar(ln)
        break;
    end

    s = strtrim(ln);
    if isempty(s)
        continue;
    end

    if s(1) ~= "#"
        break;
    end

    if ~isfinite(meta.T)
        tok = regexp(s, ['^\#\s*T\s*=\s*' num], 'tokens','once');
        if ~isempty(tok)
            meta.T = str2double(tok{1});
        end
    end

    if ~isfinite(meta.doping)
        tok = regexp(s, ['^\#\s*doping\s*=\s*' num], 'tokens','once');
        if ~isempty(tok)
            meta.doping = str2double(tok{1});
        end
    end

    if ~isfinite(meta.mu)
        tok = regexp(s, ['^\#\s*mu\s*=\s*' num], 'tokens','once');

        if isempty(tok)
            tok = regexp(s, ['^\#\s*EF\s*=\s*' num], 'tokens','once');
        end

        if ~isempty(tok)
            meta.mu = str2double(tok{1});
        end
    end
end

end

% ============================================================
% read chi table
% ============================================================
function M = read_numeric_skiphash_(fpath)

fid = fopen(fpath,"r");
if fid < 0
    M = [];
    return;
end

c = onCleanup(@() fclose(fid)); %#ok<NASGU>
rows = {};

while true

    ln = fgetl(fid);
    if ~ischar(ln)
        break;
    end

    if isempty(ln)
        continue;
    end

    if ~isempty(regexp(ln,'^\s*#','once'))
        continue;
    end

    v = sscanf(ln,'%f').';
    if isempty(v)
        continue;
    end

    rows{end+1,1} = v; %#ok<AGROW>
end

if isempty(rows)
    M = [];
    return;
end

ncol = max(cellfun(@numel, rows));
M = nan(numel(rows), ncol);

for i = 1:numel(rows)
    v = rows{i};
    M(i,1:numel(v)) = v;
end

lastFinite = find(any(isfinite(M),1), 1, "last");
if ~isempty(lastFinite)
    M = M(:,1:lastFinite);
end

end

function C = detect_cols_chi_(M)

ncol = size(M,2);

C = struct("iq",1,"jq",2,"Re",5);

if ncol >= 8
    C.iq = 2;
    C.jq = 3;
    C.Re = 6;
end

end

% ============================================================
% utilities
% ============================================================
function tf = pass_range_AND_(v, range)

if isempty(range)
    tf = true;
    return;
end

if ~isfinite(v)
    tf = false;
    return;
end

lo = min(range);
hi = max(range);
tf = (v >= lo) && (v <= hi);

end

function [mu2, dop2, chi2] = collapse_same_mu_(mu, dop, chi)

mu  = double(mu(:));
dop = double(dop(:));
chi = double(chi(:));

ok = isfinite(mu) & isfinite(dop) & isfinite(chi);
mu  = mu(ok);
dop = dop(ok);
chi = chi(ok);

if isempty(mu)
    mu2 = [];
    dop2 = [];
    chi2 = [];
    return;
end

[mu_unique, ~, ic] = unique(mu);

dop_mean = accumarray(ic, dop, [], @mean);
chi_mean = accumarray(ic, chi, [], @mean);

[mu2, ord] = sort(mu_unique);
dop2 = dop_mean(ord);
chi2 = chi_mean(ord);

end

function [x_u, y_u, z_u, w_u] = average_by_x4_(x, y, z, w)

m = isfinite(x) & isfinite(y) & isfinite(z) & isfinite(w);

x = x(m);
y = y(m);
z = z(m);
w = w(m);

if isempty(x)
    x_u = x;
    y_u = y;
    z_u = z;
    w_u = w;
    return;
end

xq = round(x, 12);
[xs, ord] = sort(xq);

ys = y(ord);
zs = z(ord);
ws = w(ord);

[ux, ~, ic] = unique(xs);

y_u = accumarray(ic, ys, [], @(v) mean(v, "omitnan"));
z_u = accumarray(ic, zs, [], @(v) mean(v, "omitnan"));
w_u = accumarray(ic, ws, [], @(v) mean(v, "omitnan"));

x_u = ux;

end

function c = polar_to_color_(p, pmin, pmax, cmap)

N = size(cmap,1);

if pmax > pmin
    t = (p - pmin) / (pmax - pmin);
else
    t = 0.5;
end

idx = 1 + round(t * (N-1));
idx = max(1, min(N, idx));

c = cmap(idx,:);

end

function cmap = blue_gradient(N)

if nargin < 1
    N = 256;
end

dark_blue  = [0.00, 0.05, 0.35];
light_blue = [0.85, 0.93, 1.00];

t = linspace(0,1,N)';
gamma = 0.6;
t = t.^gamma;

cmap = [ ...
    dark_blue(1) + (light_blue(1)-dark_blue(1)) * t, ...
    dark_blue(2) + (light_blue(2)-dark_blue(2)) * t, ...
    dark_blue(3) + (light_blue(3)-dark_blue(3)) * t ...
];

end