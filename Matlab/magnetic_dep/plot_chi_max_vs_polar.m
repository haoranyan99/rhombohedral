function H = plot_chi_max_vs_polar()

clc; close all;

% =========================
% USER SETTINGS
% =========================
iq_pick = 3;
jq_pick = -3;

T_target = 6.5;
T_tol    = 1e-8;

doping_range = [-1, 0];     % only search max chi in this doping window
use_abs_chi  = false;       % true: max |chi|, false: max chi

polar_list_meV = [];        % [] = all polar
% polar_list_meV = 0:0.5:3;

FS = 15;
LW = 2.2;
MS = 7;

default_root = "E:/rg_master/data/";

% =========================
% Select root folder
% =========================
if ~isfolder(default_root)
    default_root = pwd;
end

root_dir = uigetdir(default_root, ...
    "Select chi multipolar folder containing polar_meV*/T*/mu*/chi*.txt");

if isequal(root_dir,0)
    H = [];
    return;
end

root_dir = string(root_dir);

% =========================
% Collect chi data
% =========================
D = collect_chi_from_structure_(root_dir, iq_pick, jq_pick, T_target, T_tol);

if isempty(D)
    error("No valid chi data found.");
end

all_polar = sort(unique([D.polar_meV]));

if isempty(polar_list_meV)
    polar_use = all_polar;
else
    polar_use = polar_list_meV(:).';
end

% =========================
% Compute max chi in doping window
% =========================
polar_out = [];
chi_max_out = [];
dop_at_max_out = [];
mu_at_max_out = [];
npoint_out = [];

for ip = 1:numel(polar_use)

    p = polar_use(ip);

    idx = abs([D.polar_meV] - p) < 1e-12;
    if ~any(idx)
        fprintf("[skip] polar %.6g meV not found\n", p);
        continue;
    end

    dop = [D(idx).doping].';
    mu  = [D(idx).mu].';
    chi = [D(idx).chi].';

    if use_abs_chi
        chi_use = abs(chi);
    else
        chi_use = chi;
    end

    valid = isfinite(dop) & isfinite(mu) & isfinite(chi_use) & ...
            dop >= doping_range(1) & dop <= doping_range(2);

    if ~any(valid)
        fprintf("[skip] polar %.6g meV has no points in doping range [%.3g, %.3g]\n", ...
            p, doping_range(1), doping_range(2));
        continue;
    end

    dop_v = dop(valid);
    mu_v  = mu(valid);
    chi_v = chi_use(valid);

    [chi_max, imax] = max(chi_v);

    polar_out(end+1,1) = p; %#ok<AGROW>
    chi_max_out(end+1,1) = chi_max; %#ok<AGROW>
    dop_at_max_out(end+1,1) = dop_v(imax); %#ok<AGROW>
    mu_at_max_out(end+1,1) = mu_v(imax); %#ok<AGROW>
    npoint_out(end+1,1) = numel(chi_v); %#ok<AGROW>
end

% sort by polar
[polar_out, ord] = sort(polar_out);
chi_max_out = chi_max_out(ord);
dop_at_max_out = dop_at_max_out(ord);
mu_at_max_out = mu_at_max_out(ord);
npoint_out = npoint_out(ord);

% =========================
% Plot
% =========================
fig = figure("Color","w","Units","pixels","Position",[120 120 760 520]);
ax = axes(fig);
hold(ax,"on"); box(ax,"on");

plot(ax, polar_out, chi_max_out, "-o", ...
    "LineWidth", LW, ...
    "MarkerSize", MS, ...
    "MarkerFaceColor", "auto");

set(ax, ...
    "FontSize", FS, ...
    "LineWidth", 1.3, ...
    "TickDir", "in", ...
    "Box", "on");

grid(ax,"off");

xlabel(ax, "polar (meV)", "Interpreter","tex", "FontSize",FS);

if use_abs_chi
    ylabel(ax, "max |Re(\chi)|", "Interpreter","tex", "FontSize",FS);
else
    ylabel(ax, "max Re(\chi)", "Interpreter","tex", "FontSize",FS);
end

title(ax, sprintf("$T=%.4g$ K, $q=(%d,%d)$, doping $\\in[%.2g,%.2g]$", ...
    T_target, iq_pick, jq_pick, doping_range(1), doping_range(2)), ...
    "Interpreter","latex", "FontWeight","normal", "FontSize",FS);

% =========================
% Output
% =========================
H = struct();
H.root_dir = root_dir;
H.T_target = T_target;
H.T_tol = T_tol;
H.iq_pick = iq_pick;
H.jq_pick = jq_pick;
H.doping_range = doping_range;
H.use_abs_chi = use_abs_chi;
H.polar_meV = polar_out;
H.chi_max = chi_max_out;
H.doping_at_max = dop_at_max_out;
H.mu_at_max = mu_at_max_out;
H.npoint = npoint_out;
H.fig = fig;
H.ax = ax;

TBL = table(polar_out, chi_max_out, dop_at_max_out, mu_at_max_out, npoint_out, ...
    "VariableNames", ["polar_meV","chi_max","doping_at_max","mu_at_max","npoint"]);

disp(TBL);

end

% ============================================================
% collect real chi from:
% root / polar_meV* / T* / mu* / chi*.txt
% ============================================================
function D = collect_chi_from_structure_(root_dir, iq_pick, jq_pick, T_target, T_tol)

D = struct("polar_meV",{}, "T",{}, ...
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

            if ~isfinite(mu_val) || ~isfinite(dop_val)
                continue;
            end

            M = read_numeric_skiphash_(fpath);
            if isempty(M) || size(M,2) < 6
                continue;
            end

            C = detect_cols_chi_(M);
            id = find(M(:,C.iq) == iq_pick & M(:,C.jq) == jq_pick, 1);

            if isempty(id)
                continue;
            end

            entry.polar_meV = polar_meV;
            entry.T = Tval;
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

function C = detect_cols_chi_(M)

ncol = size(M,2);

C = struct("iq",1,"jq",2,"Re",5);

if ncol >= 8
    C.iq = 2;
    C.jq = 3;
    C.Re = 6;
end

end