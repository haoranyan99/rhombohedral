function H = plot_a1_and_chi_multipolar()

clc; close all;

% =========================
% USER SETTINGS
% =========================
iq_pick = 12;
jq_pick = -6;

T_target = 6.5;
T_tol    = 1e-8;

x_mode = "doping";          % "mu" or "doping"
use_abs_chi = false;

% 手动设定要画的真实 polar
% [] 表示画所有 polar
polar_list_meV = [];    
polar_list_meV = 0:0.5:6;

FS = 15;
LW = 2.0;

default_root = "E:/rg_master/data/";

% =========================
% Coefficients for a1
% =========================
par  = make_realchi_params(true);
coef = make_realchi_coeff(par);

% =========================
% Select root folder
% =========================
if ~isfolder(default_root)
    default_root = pwd;
end

root_dir = uigetdir(default_root, ...
    "Select chi multipolar folder containing polar_meV*/T*/mu*/chi*.txt");

if isequal(root_dir,0)
    return;
end

root_dir = string(root_dir);

% =========================
% Collect real chi data
% =========================
D = collect_chi_from_structure_(root_dir, iq_pick, jq_pick, ...
    T_target, T_tol, x_mode);

if isempty(D)
    error("No valid chi data found. Check T_target, iq/jq, x_mode, and folder structure.");
end

all_polar = sort(unique([D.polar_meV]));

if isempty(polar_list_meV)
    polar_use = all_polar;
else
    polar_use = polar_list_meV(:).';
end

fprintf("[polar used meV] ");
fprintf("%.6g ", polar_use);
fprintf("\n");

% =========================
% Color map
% =========================
pmin = min(polar_use);
pmax = max(polar_use);

cmap0 = hot(256);

% 去掉最亮的白色部分（后 20%）
cmap = cmap0(1:round(0.8*256), :);

% 再插值回256
cmap = interp1(linspace(0,1,size(cmap,1)), cmap, linspace(0,1,256));
get_color = @(p) polar_to_color_(p, pmin, pmax, cmap);

% =========================
% Prepare output
% =========================
H = struct();
H.root_dir = root_dir;
H.T_target = T_target;
H.T_tol = T_tol;
H.iq_pick = iq_pick;
H.jq_pick = jq_pick;
H.x_mode = x_mode;
H.polar_list_meV = polar_use;
H.curves = struct([]);

% =========================
% Figure 1: chi
% =========================
fig_chi = figure("Color","w","Units","pixels","Position",[120 120 900 560]);
ax_chi = axes(fig_chi);
hold(ax_chi,"on");
box(ax_chi,"on");

set(ax_chi, ...
    "FontSize",FS, ...
    "LineWidth",1.3, ...
    "TickDir","in", ...
    "Box","on");

grid(ax_chi,"off");

% =========================
% Figure 2: a1
% =========================
fig_a1 = figure("Color","w","Units","pixels","Position",[160 140 900 560]);
ax_a1 = axes(fig_a1);
hold(ax_a1,"on");
box(ax_a1,"on");

set(ax_a1, ...
    "FontSize",FS, ...
    "LineWidth",1.3, ...
    "TickDir","in", ...
    "Box","on");

grid(ax_a1,"off");

% =========================
% Loop over selected polar
% =========================
for ip = 1:numel(polar_use)

    p = polar_use(ip);

    idx = abs([D.polar_meV] - p) < 1e-12;
    if ~any(idx)
        fprintf("[skip] polar %.6g meV not found\n", p);
        continue;
    end

    x   = [D(idx).x].';
    chi = [D(idx).chi].';
    dop = [D(idx).doping].';
    mu  = [D(idx).mu].';

    if use_abs_chi
        chi_plot_raw = abs(chi);
    else
        chi_plot_raw = chi;
    end

    % average duplicate x
    [x_plot, chi_plot, dop_plot, mu_plot] = average_by_x4_(x, chi_plot_raw, dop, mu);

    % compute a1 using real chi
    a1_plot = nan(size(x_plot));
    for i = 1:numel(x_plot)
        C = coef.eval(T_target, dop_plot(i), chi_plot(i));
        a1_plot(i) = C.a1;
    end

    % sort
    [x_plot, ord] = sort(x_plot);
    chi_plot = chi_plot(ord);
    a1_plot  = a1_plot(ord);
    dop_plot = dop_plot(ord);
    mu_plot  = mu_plot(ord);

    this_color = get_color(p);

    plot(ax_chi, x_plot, chi_plot, "-", ...
        "LineWidth", LW, ...
        "Color", this_color);

    plot(ax_a1, x_plot, a1_plot, "-", ...
        "LineWidth", LW, ...
        "Color", this_color);

    H.curves(end+1).polar_meV = p; %#ok<AGROW>
    H.curves(end).x = x_plot;
    H.curves(end).chi = chi_plot;
    H.curves(end).a1 = a1_plot;
    H.curves(end).doping = dop_plot;
    H.curves(end).mu = mu_plot;
end

% =========================
% Labels
% =========================
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

title(ax_chi, sprintf("Real $\\chi$, $T=%.4g$ K, $q=(%d,%d)$", ...
    T_target, iq_pick, jq_pick), ...
    "Interpreter","latex", "FontWeight","normal", "FontSize",FS);

title(ax_a1, sprintf("$a_1$ from real $\\chi$, $T=%.4g$ K, $q=(%d,%d)$", ...
    T_target, iq_pick, jq_pick), ...
    "Interpreter","latex", "FontWeight","normal", "FontSize",FS);

yline(ax_a1, 0, "--", "LineWidth",1.1);

% =========================
% Colorbars
% =========================
colormap(ax_chi, cmap);
caxis(ax_chi, [pmin pmax]);
cb1 = colorbar(ax_chi);
cb1.Label.String = "polar (meV)";
cb1.TickDirection = "in";
cb1.LineWidth = 1.2;

colormap(ax_a1, cmap);
caxis(ax_a1, [pmin pmax]);
cb2 = colorbar(ax_a1);
cb2.Label.String = "polar (meV)";
cb2.TickDirection = "in";
cb2.LineWidth = 1.2;

H.fig_chi = fig_chi;
H.ax_chi = ax_chi;
H.fig_a1 = fig_a1;
H.ax_a1 = ax_a1;

end

% ============================================================
% collect real chi from:
% root / polar_meV* / T* / mu* / chi*.txt
% ============================================================
function D = collect_chi_from_structure_(root_dir, iq_pick, jq_pick, T_target, T_tol, x_mode)

D = struct("polar_meV",{}, "T",{}, "x",{}, ...
           "chi",{}, "mu",{}, "doping",{}, "file",{});

polar_dirs = dir(fullfile(root_dir, "polar_meV*"));
polar_dirs = polar_dirs([polar_dirs.isdir]);

fprintf("[root] %s\n", root_dir);
fprintf("[polar folders] %d\n", numel(polar_dirs));

for ip = 1:numel(polar_dirs)

    polar_name = string(polar_dirs(ip).name);
    polar_path = string(fullfile(polar_dirs(ip).folder, polar_dirs(ip).name));

    polar_meV = parse_polar_from_name_(polar_name);
    if ~isfinite(polar_meV)
        continue;
    end

    T_dirs = dir(fullfile(polar_path, "T*"));
    T_dirs = T_dirs([T_dirs.isdir]);

    for it = 1:numel(T_dirs)

        T_name = string(T_dirs(it).name);
        T_path = string(fullfile(T_dirs(it).folder, T_dirs(it).name));

        Tval = parse_T_from_name_(T_name);
        if ~isfinite(Tval)
            continue;
        end

        if abs(Tval - T_target) > T_tol
            continue;
        end

        chi_files = dir(fullfile(T_path, "**", "chi*.txt"));

        for k = 1:numel(chi_files)

            fpath = string(fullfile(chi_files(k).folder, chi_files(k).name));

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
            if isempty(M) || size(M,2) < 6
                continue;
            end

            C = detect_chi_columns(M, fpath, iq_pick, jq_pick);
            id = find(M(:,C.iq) == iq_pick & M(:,C.jq) == jq_pick, 1);

            if isempty(id)
                continue;
            end

            entry.polar_meV = polar_meV;
            entry.T = Tval;
            entry.x = xv;
            entry.chi = M(id, C.Re);
            entry.mu = mu_val;
            entry.doping = dop_val;
            entry.file = fpath;

            D(end+1) = entry; %#ok<AGROW>
        end
    end
end

fprintf("[valid points] %d\n", numel(D));

end

% ============================================================
% parsers
% ============================================================
function polar_meV = parse_polar_from_name_(name)

polar_meV = NaN;
tok = regexp(char(name), ...
    '^polar_meV([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)$', ...
    'tokens', 'once');

if ~isempty(tok)
    polar_meV = str2double(tok{1});
end

end

function Tval = parse_T_from_name_(name)

Tval = NaN;
tok = regexp(char(name), ...
    '^T([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)$', ...
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

meta = struct("doping",NaN,"mu",NaN);

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

% ============================================================
% average duplicate x
% ============================================================
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

% ============================================================
% color helper
% ============================================================
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
