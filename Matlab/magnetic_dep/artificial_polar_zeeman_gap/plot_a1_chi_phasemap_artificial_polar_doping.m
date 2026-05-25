function H = plot_a1_chi_phasemap_artificial_polar_doping(root_dir)
% ============================================================
% USER SETTINGS
% ============================================================
T_target = 6.4;
T_tol    = 1e-8;

iq_pick = 12;
jq_pick = -6;

polar_list_meV = 0:0.1:2.5;

% Manual display range. Use [] to auto-select from the common valid range.
doping_range_plot = [-1.9, 0];

% Optional read filter before building artificial-polar averages.
doping_range_read = [];
mu_range = [];
doping_N = 350;

FS = 15;
chi_white_pos = 0.65;    % smaller value gives more white-to-red levels
a1_white_pos = 0.40;     % larger value gives more red-to-white levels

% Re(chi) flipped-hot color settings.
chi_hot_white_value = 0.02;      % Re(chi) below this stays near white
chi_hot_white_color_pos = 0.06;  % position in flip-hot colormap for white_frac
chi_hot_red_value = 0.0695;       % Re(chi) value shown as red
chi_hot_red_color_pos = 0.65;    % red position in flip-hot colormap
chi_hot_nticks = 4;
chi_cap_doping_window = [-1.6, -1.4];
chi_cap_value = 0.069;

save_outputs = true;

% ============================================================
% Locate chi root
% ============================================================
if nargin < 1 || strlength(string(root_dir)) == 0
    cand = [
        string(fullfile(pwd, "data", "chi_12_m6"))
        string(fullfile(pwd, "data", "chi_12m6"))
        "E:/rg_master/rhombohedral/data/chi_12_m6"
        "E:/rg_master/rhombohedral/data/chi_12m6"
    ];

    root_dir = "";
    for i = 1:numel(cand)
        if isfolder(cand(i))
            root_dir = cand(i);
            break;
        end
    end

    if strlength(root_dir) == 0
        pick = uigetdir(pwd, "Select chi_12_m6 folder containing T*/mu*/chi*.txt");
        if isequal(pick, 0)
            return;
        end
        root_dir = string(pick);
    end
else
    root_dir = string(root_dir);
end

if ~isfolder(root_dir)
    error("Chi root folder not found: %s", root_dir);
end

% ============================================================
% Coefficients for a1
% ============================================================
par  = make_realchi_params(true);
coef = make_realchi_coeff(par);

% ============================================================
% Collect baseline chi curve
% ============================================================
[mu0, dop0, chi0] = collect_curve_from_folder_( ...
    root_dir, T_target, T_tol, iq_pick, jq_pick, mu_range, doping_range_read);

if numel(mu0) < 3
    error("Baseline curve has too few valid points.");
end

mu0  = double(mu0(:));
dop0 = double(dop0(:));
chi0 = double(chi0(:));

[mu0, ord] = sort(mu0);
dop0 = dop0(ord);
chi0 = chi0(ord);

[mu0, dop0, chi0] = collapse_same_mu_(mu0, dop0, chi0);

Fchi = griddedInterpolant(mu0, chi0, "linear", "none");

polar_use = unique(double(polar_list_meV(:)));
polar_use = polar_use(isfinite(polar_use));
polar_use = sort(polar_use);

curves = struct("polar_meV",{}, "dop",{}, "chi",{}, "a1",{});

% ============================================================
% Build artificial-polar curves
% ============================================================
for ip = 1:numel(polar_use)
    p = polar_use(ip);
    dmu = 0.5 * p * 1e-3;  % meV -> eV, half splitting

    chi_plus  = Fchi(mu0 + dmu);
    chi_minus = Fchi(mu0 - dmu);
    chi_art   = 0.5 * (chi_plus + chi_minus);

    ok = isfinite(dop0) & isfinite(chi_art);
    if nnz(ok) < 3
        fprintf("[skip] artificial polar %.6g meV has too few valid points\n", p);
        continue;
    end

    dop = dop0(ok);
    chi = chi_art(ok);

    [dop, chi] = average_by_x2_(dop, chi);
    [dop, ord] = sort(dop);
    chi = chi(ord);

    a1 = nan(size(chi));
    for i = 1:numel(chi)
        C = coef.eval(T_target, dop(i), chi(i));
        a1(i) = C.a1;
    end

    curves(end+1).polar_meV = p; %#ok<AGROW>
    curves(end).dop = dop;
    curves(end).chi = chi;
    curves(end).a1 = a1;
end

if isempty(curves)
    error("No valid artificial-polar curves were built.");
end

% ============================================================
% Common valid doping grid
% ============================================================
dop_common_min = -inf;
dop_common_max = inf;
for i = 1:numel(curves)
    dop_common_min = max(dop_common_min, min(curves(i).dop));
    dop_common_max = min(dop_common_max, max(curves(i).dop));
end

if isempty(doping_range_plot)
    dop_min = dop_common_min;
    dop_max = dop_common_max;
else
    dop_min = max(min(doping_range_plot), dop_common_min);
    dop_max = min(max(doping_range_plot), dop_common_max);
end

if ~(isfinite(dop_min) && isfinite(dop_max) && dop_max > dop_min)
    error("Cannot build common doping range. Common valid range is [%.6g, %.6g].", ...
        dop_common_min, dop_common_max);
end

dop_grid = linspace(dop_min, dop_max, doping_N);
polar_grid = [curves.polar_meV].';

Z_chi = nan(numel(curves), doping_N);
Z_a1  = nan(numel(curves), doping_N);

for i = 1:numel(curves)
    [xd, ia] = unique(curves(i).dop(:));
    chi = curves(i).chi(:);
    a1  = curves(i).a1(:);

    Z_chi(i,:) = interp1(xd, chi(ia), dop_grid, "linear", NaN);
    Z_a1(i,:)  = interp1(xd, a1(ia),  dop_grid, "linear", NaN);
end

% ============================================================
% Plot
% ============================================================
Z_chi_plot = cap_value_in_doping_window_( ...
    Z_chi, dop_grid, chi_cap_doping_window, chi_cap_value);

fig_chi = plot_chi_hot_map_(dop_grid, polar_grid, Z_chi_plot, ...
    "Re(\chi)", ...
    sprintf("Artificial-polar phase map of Re(\\chi), T=%.4g K, q=(%d,%d)", ...
        T_target, iq_pick, jq_pick), ...
    FS, chi_hot_white_value, chi_hot_white_color_pos, ...
    chi_hot_red_value, chi_hot_red_color_pos, chi_hot_nticks);

fig_a1 = plot_one_map_(dop_grid, polar_grid, Z_a1, ...
    "a_1", ...
    sprintf("Artificial-polar phase map of a_1, T=%.4g K, q=(%d,%d)", ...
        T_target, iq_pick, jq_pick), ...
    FS, 0, true, a1_white_pos, true);

H = struct();
H.root_dir = root_dir;
H.T_target = T_target;
H.T_tol = T_tol;
H.iq_pick = iq_pick;
H.jq_pick = jq_pick;
H.polar_meV = polar_grid;
H.doping = dop_grid(:);
H.doping_common_range = [dop_common_min, dop_common_max];
H.doping_range_plot = [dop_min, dop_max];
H.Z_chi = Z_chi;
H.Z_chi_plot = Z_chi_plot;
H.Z_a1 = Z_a1;
H.curves = curves;
H.fig_chi = fig_chi;
H.fig_a1 = fig_a1;

if save_outputs
    out_dir = fullfile(root_dir, "plot");
    if ~exist(out_dir, "dir")
        mkdir(out_dir);
    end

    tag = sprintf("T%.4g_q%d_%d", T_target, iq_pick, jq_pick);
    png_chi = fullfile(out_dir, "phasemap_chi_artificialPolar_doping_" + tag + ".png");
    png_a1  = fullfile(out_dir, "phasemap_a1_artificialPolar_doping_"  + tag + ".png");
    mat_out = fullfile(out_dir, "phasemap_a1_chi_artificialPolar_doping_" + tag + ".mat");

    exportgraphics(fig_chi, png_chi, "Resolution", 300);
    exportgraphics(fig_a1,  png_a1,  "Resolution", 300);
    Hsave = rmfield(H, ["fig_chi","fig_a1"]);
    save(mat_out, "Hsave");

    H.png_chi = string(png_chi);
    H.png_a1 = string(png_a1);
    H.mat = string(mat_out);

    fprintf("Saved chi map: %s\n", png_chi);
    fprintf("Saved a1  map: %s\n", png_a1);
    fprintf("Saved data  : %s\n", mat_out);
end

fprintf("Root: %s\n", root_dir);
fprintf("Valid polar curves: %d\n", numel(curves));
fprintf("Common valid doping range: [%.6g, %.6g]\n", dop_common_min, dop_common_max);
fprintf("Doping grid: [%.6g, %.6g], N=%d\n", dop_grid(1), dop_grid(end), numel(dop_grid));

end

% ============================================================
% Plot helper
% ============================================================
function fig = plot_chi_hot_map_(dop_grid, polar_grid, Z, cb_label, ttl, FS, ...
    white_value, white_color_pos, red_value, red_color_pos, nticks)

fig = figure("Color","w","Units","pixels","Position",[120 120 880 640]);
ax = axes(fig);

zmin = 0;
zmax = max(Z(:), [], "omitnan");

if ~(isfinite(zmax) && zmax > zmin)
    zmax = zmin + eps;
end

if nargin < 7 || ~isfinite(white_value)
    white_value = 0.02;
end
if nargin < 8 || ~isfinite(white_color_pos)
    white_color_pos = 0.06;
end
if nargin < 9 || ~isfinite(red_value)
    red_value = 0.069;
end
if nargin < 10 || ~isfinite(red_color_pos)
    red_color_pos = 0.65;
end
if nargin < 11 || ~isfinite(nticks)
    nticks = 4;
end

red_value = min(max(red_value, zmin), zmax);
red_frac = (red_value - zmin) / (zmax - zmin);
white_value = min(max(white_value, zmin), zmax);
white_frac = (white_value - zmin) / (zmax - zmin);

imagesc(ax, dop_grid, polar_grid, Z);
set(ax, "YDir","normal");

colormap(ax, flipped_hot_with_red_anchor_( ...
    256, white_frac, white_color_pos, red_frac, red_color_pos));
clim(ax, [zmin, zmax]);

cb = colorbar(ax);
cb.Label.String = cb_label;
cb.Label.Interpreter = "tex";
cb.TickDirection = "in";
cb.LineWidth = 1.2;

ticks = linspace(zmin, zmax, max(2, round(nticks)));
cb.Ticks = ticks;
cb.TickLabels = cellstr(compose("%.1g", ticks));

fprintf("%s hot color range = [%.12g, %.12g], red at %.12g\n", ...
    cb_label, zmin, zmax, red_value);

xlabel(ax, "doping (10^{12} cm^{-2})", "Interpreter","tex");
ylabel(ax, "artificial polar (meV)", "Interpreter","tex");
title(ax, ttl, "Interpreter","tex", "FontWeight","normal");

set(ax, ...
    "FontSize", FS, ...
    "LineWidth", 1.3, ...
    "TickDir", "in", ...
    "Box", "on");

grid(ax, "off");

end

function Zout = cap_value_in_doping_window_(Z, dop_grid, dop_window, cap_value)

Zout = Z;

if isempty(dop_window) || numel(dop_window) < 2 || ~isfinite(cap_value)
    return;
end

dop_lo = min(dop_window);
dop_hi = max(dop_window);
mask = dop_grid >= dop_lo & dop_grid <= dop_hi;

if any(mask)
    Zout(:, mask) = min(Zout(:, mask), cap_value);
end

end

function fig = plot_one_map_(dop_grid, polar_grid, Z, cb_label, ttl, FS, white_value, flip_colors, white_pos_manual, integer_endpoints)

fig = figure("Color","w","Units","pixels","Position",[120 120 880 640]);
ax = axes(fig);

zmin = min(Z(:), [], "omitnan");
zmax = max(Z(:), [], "omitnan");

if nargin < 7 || ~isfinite(white_value)
    white_value = 0;
end
if nargin < 8
    flip_colors = false;
end
if nargin < 9
    white_pos_manual = [];
end
if nargin < 10
    integer_endpoints = false;
end

if ~(isfinite(zmin) && isfinite(zmax))
    zmin = -1;
    zmax = 1;
end

zmin = min(zmin, white_value);
zmax = max(zmax, white_value);

if integer_endpoints
    zmin = floor(zmin);
    zmax = ceil(zmax);
end

if ~(zmax > zmin)
    pad = max(1e-12, abs(zmin) * 1e-6);
    zmin = zmin - pad;
    zmax = zmax + pad;
end

if ~isempty(white_pos_manual) && isfinite(white_pos_manual)
    white_pos = white_pos_manual;
else
    white_pos = (white_value - zmin) / (zmax - zmin);
end
white_pos = min(max(white_pos, 0), 1);

if ~isempty(white_pos_manual) && isfinite(white_pos_manual)
    Zplot = nonlinear_color_value_(Z, zmin, white_value, zmax, white_pos);
    imagesc(ax, dop_grid, polar_grid, Zplot);
    color_limits = [0, 1];
    color_ticks = [0, white_pos, 1];
    color_tick_labels = compose("%.4g", [zmin, white_value, zmax]);
else
    imagesc(ax, dop_grid, polar_grid, Z);
    color_limits = [zmin, zmax];
    color_ticks = unique([zmin, white_value, zmax]);
    color_tick_labels = compose("%.4g", color_ticks);
end
set(ax, "YDir","normal");

cmap = red_white_blue_skewed(256, white_pos, flip_colors);
colormap(ax, cmap);
clim(ax, color_limits);

cb = colorbar(ax);
cb.Label.String = cb_label;
cb.Label.Interpreter = "tex";
cb.TickDirection = "in";
cb.LineWidth = 1.2;
cb.Ticks = color_ticks;
cb.TickLabels = color_tick_labels;

fprintf("%s color range = [%.12g, %.12g], white at %.12g\n", ...
    cb_label, zmin, zmax, white_value);

xlabel(ax, "doping (10^{12} cm^{-2})", "Interpreter","tex");
ylabel(ax, "artificial polar (meV)", "Interpreter","tex");
title(ax, ttl, "Interpreter","tex", "FontWeight","normal");

set(ax, ...
    "FontSize", FS, ...
    "LineWidth", 1.3, ...
    "TickDir", "in", ...
    "Box", "on");

grid(ax, "off");

end

function Zplot = nonlinear_color_value_(Z, zmin, white_value, zmax, white_pos)

Zplot = nan(size(Z));

lo = isfinite(Z) & Z <= white_value;
hi = isfinite(Z) & Z > white_value;

den_lo = white_value - zmin;
den_hi = zmax - white_value;

if den_lo > 0
    Zplot(lo) = white_pos * (Z(lo) - zmin) / den_lo;
else
    Zplot(lo) = white_pos;
end

if den_hi > 0
    Zplot(hi) = white_pos + (1 - white_pos) * (Z(hi) - white_value) / den_hi;
else
    Zplot(hi) = white_pos;
end

end

function cmap = red_white_blue_skewed(N, white_pos, flip_colors)

if nargin < 1
    N = 256;
end
if nargin < 2
    white_pos = 0.5;
end
if nargin < 3
    flip_colors = false;
end

blue  = [45, 75, 145] / 255;
white = [1, 1, 1];
red   = [175, 20, 20] / 255;

if flip_colors
    low_color = red;
    high_color = blue;
else
    low_color = blue;
    high_color = red;
end

n_blue = max(2, round(N * white_pos));
n_red  = N - n_blue;

c1 = [linspace(low_color(1), white(1), n_blue)', ...
      linspace(low_color(2), white(2), n_blue)', ...
      linspace(low_color(3), white(3), n_blue)'];

c2 = [linspace(white(1), high_color(1), n_red)', ...
      linspace(white(2), high_color(2), n_red)', ...
      linspace(white(3), high_color(3), n_red)'];

cmap = [c1; c2];

end

function cmap = flipped_hot_with_red_anchor_(N, white_value_frac, white_color_pos, red_value_frac, red_color_pos)

if nargin < 1
    N = 256;
end
if nargin < 2
    white_value_frac = 0.30;
end
if nargin < 3
    white_color_pos = 0.06;
end
if nargin < 4
    red_value_frac = 0.9;
end
if nargin < 5
    red_color_pos = 0.65;
end

white_value_frac = min(max(white_value_frac, 1e-6), 1 - 2e-6);
white_color_pos = min(max(white_color_pos, 1e-6), 1 - 2e-6);
red_value_frac = min(max(red_value_frac, 1e-6), 1 - 1e-6);
if red_value_frac <= white_value_frac
    red_value_frac = min(1 - 1e-6, white_value_frac + 1e-6);
end
red_color_pos = min(max(red_color_pos, white_color_pos + 1e-6), 1 - 1e-6);

base = flipud(hot(max(1024, 4 * N)));
r = linspace(0, 1, N).';
u = nan(size(r));

lo = r <= white_value_frac;
mid = r > white_value_frac & r <= red_value_frac;
hi = r > red_value_frac;
u(lo) = white_color_pos * r(lo) / white_value_frac;
u(mid) = white_color_pos + (red_color_pos - white_color_pos) * ...
    (r(mid) - white_value_frac) / (red_value_frac - white_value_frac);
u(hi) = red_color_pos + (1 - red_color_pos) * ...
    (r(hi) - red_value_frac) / (1 - red_value_frac);

xbase = linspace(0, 1, size(base, 1)).';
cmap = interp1(xbase, base, u, "linear");

end

% ============================================================
% collect baseline chi curve
% ============================================================
function [mu0, dop0, chi0] = collect_curve_from_folder_(folder, T_target, T_tol, ...
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

    C = detect_chi_columns(M, fpath, iq_pick, jq_pick);
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
% parsers / utilities
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
        tok = regexp(s, ['^\#\s*T(?:_K)?\s*=\s*' num], 'tokens','once');
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
        tok = regexp(s, ['^\#\s*(?:mu|EF|E_F)\s*=\s*' num], 'tokens','once');
        if ~isempty(tok)
            meta.mu = str2double(tok{1});
        end
    end
end

end

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

function [x_u, y_u] = average_by_x2_(x, y)

m = isfinite(x) & isfinite(y);
x = x(m);
y = y(m);

if isempty(x)
    x_u = x;
    y_u = y;
    return;
end

xq = round(x, 12);
[xs, ord] = sort(xq);
ys = y(ord);

[x_u, ~, ic] = unique(xs);
y_u = accumarray(ic, ys, [], @(v) mean(v, "omitnan"));

end
