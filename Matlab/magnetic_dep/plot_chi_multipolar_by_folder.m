function plot_chi_multipolar_by_folder()
clc; close all;

%% user settings
iq_pick = 3;
jq_pick = -3;

T_target = 6.5;
T_tol    = 1e-8;

x_mode = "mu";       % "mu" or "doping"
use_abs_chi = false;

FS = 15;

%% select root
root_dir = uigetdir(pwd, "Select chi_400_3_-3multipolar folder");
if isequal(root_dir,0), return; end
root_dir = string(root_dir);

D = collect_chi_from_structure_(root_dir, iq_pick, jq_pick, T_target, T_tol, x_mode);

if isempty(D)
    error("No valid chi data found. Check T_target, iq/jq, x_mode, and folder structure.");
end

polar_list = sort(unique([D.polar_meV]));

%% plot
figure("Color","w","Units","pixels","Position",[120 120 900 560]);
ax = axes; hold(ax,"on"); box(ax,"on"); grid(ax,"on");
set(ax,"FontSize",FS,"LineWidth",1.2,"TickDir","out");

cmap = hot(max(numel(polar_list)+2,4));
cmap = cmap(1:numel(polar_list),:);

for ip = 1:numel(polar_list)
    p = polar_list(ip);
    idx = [D.polar_meV] == p;

    x = [D(idx).x].';
    y = [D(idx).chi].';

    if use_abs_chi
        y = abs(y);
    end

    [x, y] = average_by_x_(x, y);

    plot(ax, x, y, "o-", ...
        "LineWidth",1.8, ...
        "MarkerSize",4.5, ...
        "Color",cmap(ip,:), ...
        "DisplayName",sprintf("polar = %.3f meV", p));
end

if x_mode == "mu"
    xlabel(ax,"\mu (eV)","Interpreter","tex","FontSize",FS);
else
    xlabel(ax,"doping (10^{12} cm^{-2})","Interpreter","tex","FontSize",FS);
end

if use_abs_chi
    ylabel(ax,"|Re(\chi)|","Interpreter","tex","FontSize",FS);
else
    ylabel(ax,"Re(\chi)","Interpreter","tex","FontSize",FS);
end

title(ax, sprintf("T = %.4g K, q = (%d,%d)", T_target, iq_pick, jq_pick), ...
    "Interpreter","none", ...
    "FontWeight","normal", ...
    "FontSize",FS);

legend(ax,"Location","best","Interpreter","none");

end

%% ============================================================
% collect according to exact structure:
% root / polar_meV* / T* / mu* / chi*.txt
%% ============================================================
function D = collect_chi_from_structure_(root_dir, iq_pick, jq_pick, T_target, T_tol, x_mode)

D = struct("polar_meV",{}, "T",{}, "x",{}, "chi",{}, "file",{});

polar_dirs = dir(fullfile(root_dir, "polar_meV*"));
polar_dirs = polar_dirs([polar_dirs.isdir]);

fprintf("[root] %s\n", root_dir);
fprintf("[polar folders] %d\n", numel(polar_dirs));

for ip = 1:numel(polar_dirs)
    polar_name = string(polar_dirs(ip).name);
    polar_path = string(fullfile(polar_dirs(ip).folder, polar_dirs(ip).name));

    polar_meV = parse_polar_from_name_(polar_name);
    if ~isfinite(polar_meV)
        fprintf("[skip polar] cannot parse: %s\n", polar_name);
        continue;
    end

    T_dirs = dir(fullfile(polar_path, "T*"));
    T_dirs = T_dirs([T_dirs.isdir]);

    for it = 1:numel(T_dirs)
        T_name = string(T_dirs(it).name);
        T_path = string(fullfile(T_dirs(it).folder, T_dirs(it).name));

        Tval = parse_T_from_name_(T_name);
        if ~isfinite(Tval)
            fprintf("[skip T] cannot parse: %s\n", T_name);
            continue;
        end

        if abs(Tval - T_target) > T_tol
            continue;
        end

        chi_files = dir(fullfile(T_path, "**", "chi*.txt"));

        for k = 1:numel(chi_files)
            fpath = string(fullfile(chi_files(k).folder, chi_files(k).name));

            H = parse_header_mu_doping_(fpath);

            if x_mode == "mu"
                xv = H.mu;
                if ~isfinite(xv)
                    xv = parse_mu_from_path_(fpath);
                end
            else
                xv = H.doping;
            end

            if ~isfinite(xv)
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
            entry.x = xv;
            entry.chi = M(id, C.Re);
            entry.file = fpath;

            D(end+1) = entry; %#ok<AGROW>
        end
    end
end

fprintf("[valid points] %d\n", numel(D));

end

%% ============================================================
% parsers
%% ============================================================
function polar_meV = parse_polar_from_name_(name)
polar_meV = NaN;
tok = regexp(char(name), '^polar_meV([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)$', 'tokens', 'once');
if ~isempty(tok)
    polar_meV = str2double(tok{1});
end
end

function Tval = parse_T_from_name_(name)
Tval = NaN;
tok = regexp(char(name), '^T([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)$', 'tokens', 'once');
if ~isempty(tok)
    Tval = str2double(tok{1});
end
end

function mu = parse_mu_from_path_(fpath)
mu = NaN;
tok = regexp(char(fpath), '[\\/]+mu([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)[\\/]+', 'tokens', 'once');
if ~isempty(tok)
    mu = str2double(tok{1});
end
end

function meta = parse_header_mu_doping_(fpath)
meta = struct("doping",NaN,"mu",NaN);

fid = fopen(fpath,"r");
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';

for t = 1:4000
    ln = fgetl(fid);
    if ~ischar(ln), break; end

    s = strtrim(ln);
    if isempty(s), continue; end
    if s(1) ~= "#", break; end

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

%% ============================================================
% read chi table
%% ============================================================
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
    if ~ischar(ln), break; end
    if isempty(ln), continue; end
    if ~isempty(regexp(ln,'^\s*#','once')), continue; end

    v = sscanf(ln,'%f').';
    if isempty(v), continue; end

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

function [x_u, y_u] = average_by_x_(x, y)
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

[ux, ~, ic] = unique(xs);
y_u = accumarray(ic, ys, [], @(v) mean(v, "omitnan"));
x_u = ux;
end