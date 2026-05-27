function H = plot_chi_branch_vs_mu()
clc; close all;

% =========================
% USER SETTINGS
% =========================
iq_pick = 4;
jq_pick = 0;

T_target = 6.4;
T_tol = 1e-8;

x_mode = "doping";  % "mu" or "doping"
use_abs_chi = false;

% "average", "valley_plus_spin_up", "valley_plus_spin_down",
% "valley_minus_spin_up", or "valley_minus_spin_down"
branch_choice = "average";

% [] means use all B folders found under the selected chi root.
B_list_T = [];

FS = 15;
LW = 2.0;

default_root = "E:/rg_master/rhombohedral/data/";

% =========================
% Select root folder
% =========================
if ~isfolder(default_root)
    default_root = pwd;
end

root_dir = uigetdir(default_root, ...
    "Select chi run root, e.g. data/chi_YYYY-MM-DD_HHMMSS");

if isequal(root_dir, 0)
    H = [];
    return;
end

root_dir = string(root_dir);

% =========================
% Collect chi data
% =========================
D = collect_chi_from_B_tree_(root_dir, branch_choice, iq_pick, jq_pick, ...
    T_target, T_tol, x_mode);

if isempty(D)
    error("No valid chi data found for branch '%s'.", branch_choice);
end

all_B = sort(unique([D.B_T]));

if isempty(B_list_T)
    B_use = all_B;
else
    B_use = B_list_T(:).';
end

fprintf("[branch] %s\n", branch_choice);
fprintf("[B used T] ");
fprintf("%.6g ", B_use);
fprintf("\n");

cmap = hot(256);
Bmin = min(B_use);
Bmax = max(B_use);
get_color = @(B) B_to_color_(B, Bmin, Bmax, cmap);

% =========================
% Plot
% =========================
fig = figure("Color", "w", "Units", "pixels", "Position", [120 120 860 540]);
ax = axes(fig);
hold(ax, "on");
box(ax, "on");

H = struct();
H.root_dir = root_dir;
H.branch_choice = branch_choice;
H.T_target = T_target;
H.T_tol = T_tol;
H.iq_pick = iq_pick;
H.jq_pick = jq_pick;
H.x_mode = x_mode;
H.B_list_T = B_use;
H.curves = struct([]);

for iB = 1:numel(B_use)
    B = B_use(iB);
    idx = abs([D.B_T] - B) < 1e-12;
    if ~any(idx)
        fprintf("[skip] B %.6g T not found\n", B);
        continue;
    end

    x = [D(idx).x].';
    chi = [D(idx).chi].';
    dop = [D(idx).doping].';
    mu = [D(idx).mu].';

    if use_abs_chi
        chi_plot_raw = abs(chi);
    else
        chi_plot_raw = chi;
    end

    [x_plot, chi_plot, dop_plot, mu_plot] = average_by_x4_(x, chi_plot_raw, dop, mu);
    [x_plot, ord] = sort(x_plot);
    chi_plot = chi_plot(ord);
    dop_plot = dop_plot(ord);
    mu_plot = mu_plot(ord);

    plot(ax, x_plot, chi_plot, "-", ...
        "LineWidth", LW, ...
        "Color", get_color(B), ...
        "HandleVisibility", "off");

    H.curves(end + 1).B_T = B; %#ok<AGROW>
    H.curves(end).x = x_plot;
    H.curves(end).chi = chi_plot;
    H.curves(end).doping = dop_plot;
    H.curves(end).mu = mu_plot;
end

set(ax, ...
    "FontSize", FS, ...
    "LineWidth", 1.3, ...
    "TickDir", "in", ...
    "Box", "on");

grid(ax, "off");

if x_mode == "mu"
    xlab = "\mu (eV)";
else
    xlab = "doping (10^{12} cm^{-2})";
end

xlabel(ax, xlab, "Interpreter", "tex", "FontSize", FS);

if use_abs_chi
    ylabel(ax, "|Re(\chi)|", "Interpreter", "tex", "FontSize", FS);
else
    ylabel(ax, "Re(\chi)", "Interpreter", "tex", "FontSize", FS);
end

title(ax, sprintf("%s, $T=%.4g$ K, $q=(%d,%d)$", ...
    branch_choice, T_target, iq_pick, jq_pick), ...
    "Interpreter", "latex", "FontWeight", "normal", "FontSize", FS);

colormap(ax, cmap);
caxis(ax, [Bmin Bmax]);
cb = colorbar(ax);
cb.Label.String = "B [T]";
cb.TickDirection = "in";
cb.LineWidth = 1.2;

H.fig = fig;
H.ax = ax;

end

% ============================================================
% Expected tree:
% root / D* / B* / T* / mu* / chi_<branch>_*.txt
% root may also be a D* or B* subfolder.
% ============================================================
function D = collect_chi_from_B_tree_(root_dir, branch_choice, iq_pick, jq_pick, ...
    T_target, T_tol, x_mode)

D = struct("B_T", {}, "T", {}, "x", {}, ...
           "chi", {}, "mu", {}, "doping", {}, "branch", {}, "file", {});

pattern = chi_pattern_from_branch_(branch_choice);
chi_files = dir(fullfile(root_dir, "**", pattern));

fprintf("[root] %s\n", root_dir);
fprintf("[pattern] %s\n", pattern);
fprintf("[candidate files] %d\n", numel(chi_files));

for k = 1:numel(chi_files)
    fpath = string(fullfile(chi_files(k).folder, chi_files(k).name));

    Tval = parse_T_from_path_(fpath);
    if ~isfinite(Tval) || abs(Tval - T_target) > T_tol
        continue;
    end

    Bval = parse_B_from_path_(fpath);
    if ~isfinite(Bval)
        continue;
    end

    H = parse_header_mu_doping_(fpath);

    mu_val = H.mu;
    if ~isfinite(mu_val)
        mu_val = parse_mu_from_path_(fpath);
    end

    dop_val = H.doping;

    if x_mode == "mu"
        xv = mu_val;
    else
        xv = dop_val;
    end

    if ~isfinite(xv) || ~isfinite(mu_val) || ~isfinite(dop_val)
        continue;
    end

    M = read_numeric_skiphash_(fpath);
    if isempty(M) || size(M, 2) < 6
        continue;
    end

    C = detect_chi_columns(M, fpath, iq_pick, jq_pick);
    id = find(M(:, C.iq) == iq_pick & M(:, C.jq) == jq_pick, 1);

    if isempty(id)
        continue;
    end

    entry.B_T = Bval;
    entry.T = Tval;
    entry.x = xv;
    entry.chi = M(id, C.Re);
    entry.mu = mu_val;
    entry.doping = dop_val;
    entry.branch = parse_branch_from_name_(string(chi_files(k).name));
    entry.file = fpath;

    D(end + 1) = entry; %#ok<AGROW>
end

if branch_choice == "average" && ~isempty(D)
    D = average_branch_entries_(D);
end

fprintf("[valid points] %d\n", numel(D));
end

function pattern = chi_pattern_from_branch_(branch_choice)
branch_choice = string(branch_choice);
if branch_choice == "average"
    pattern = "chi_valley_*_spin_*_*.txt";
else
    pattern = "chi_" + branch_choice + "_*.txt";
end
pattern = char(pattern);
end

function branch = parse_branch_from_name_(name)
branch = "";
tok = regexp(char(name), ...
    '^chi_(valley_(?:plus|minus)_spin_(?:up|down))_', ...
    'tokens', 'once');
if ~isempty(tok)
    branch = string(tok{1});
end
end

function D_avg = average_branch_entries_(D)
keys = strings(1, numel(D));
for i = 1:numel(D)
    keys(i) = sprintf("B%.12g_T%.12g_mu%.12g", ...
        D(i).B_T, D(i).T, D(i).mu);
end

[ukeys, ~, ic] = unique(keys);
D_avg = struct("B_T", {}, "T", {}, "x", {}, ...
               "chi", {}, "mu", {}, "doping", {}, "branch", {}, "file", {});

for i = 1:numel(ukeys)
    idx = find(ic == i);
    branches = unique([D(idx).branch]);

    if numel(branches) < 4
        fprintf("[warn] averaging %d branch file(s) for %s\n", ...
            numel(branches), ukeys(i));
    end

    entry = D(idx(1));
    entry.x = mean([D(idx).x], "omitnan");
    entry.chi = mean([D(idx).chi], "omitnan");
    entry.doping = mean([D(idx).doping], "omitnan");
    entry.mu = mean([D(idx).mu], "omitnan");
    entry.branch = "average";
    entry.file = strjoin([D(idx).file], ";");

    D_avg(end + 1) = entry; %#ok<AGROW>
end
end

function B = parse_B_from_path_(fpath)
B = NaN;
tok = regexp(char(fpath), ...
    '[\\/]+B([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)T[\\/]+', ...
    'tokens', 'once');
if ~isempty(tok)
    B = str2double(tok{1});
end
end

function Tval = parse_T_from_path_(fpath)
Tval = NaN;
tok = regexp(char(fpath), ...
    '[\\/]+T([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)[\\/]+', ...
    'tokens', 'once');
if ~isempty(tok)
    Tval = str2double(tok{1});
end
end

function mu = parse_mu_from_path_(fpath)
mu = NaN;
tok = regexp(char(fpath), ...
    '[\\/]+mu([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)[\\/]+', ...
    'tokens', 'once');
if ~isempty(tok)
    mu = str2double(tok{1});
end
end

function meta = parse_header_mu_doping_(fpath)
meta = struct("doping", NaN, "mu", NaN);

fid = fopen(fpath, "r");
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

    if ~isfinite(meta.doping)
        tok = regexp(s, ['^\#\s*doping\s*=\s*' num], 'tokens', 'once');
        if ~isempty(tok)
            meta.doping = str2double(tok{1});
        end
    end

    if ~isfinite(meta.mu)
        tok = regexp(s, ['^\#\s*mu\s*=\s*' num], 'tokens', 'once');
        if isempty(tok)
            tok = regexp(s, ['^\#\s*EF\s*=\s*' num], 'tokens', 'once');
        end
        if ~isempty(tok)
            meta.mu = str2double(tok{1});
        end
    end
end
end

function M = read_numeric_skiphash_(fpath)
fid = fopen(fpath, "r");
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

    if isempty(ln) || ~isempty(regexp(ln, '^\s*#', 'once'))
        continue;
    end

    v = sscanf(ln, '%f').';
    if isempty(v)
        continue;
    end

    rows{end + 1, 1} = v; %#ok<AGROW>
end

if isempty(rows)
    M = [];
    return;
end

ncol = max(cellfun(@numel, rows));
M = nan(numel(rows), ncol);

for i = 1:numel(rows)
    v = rows{i};
    M(i, 1:numel(v)) = v;
end

lastFinite = find(any(isfinite(M), 1), 1, "last");
if ~isempty(lastFinite)
    M = M(:, 1:lastFinite);
end
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

function c = B_to_color_(B, Bmin, Bmax, cmap)
N = size(cmap, 1);

if Bmax > Bmin
    t = (B - Bmin) / (Bmax - Bmin);
else
    t = 0.5;
end

idx = 1 + round(t * (N - 1));
idx = max(1, min(N, idx));
c = cmap(idx, :);
end
