function plot_chi()


% =========================
% User settings
% =========================
default_root = "E:/rg_master/rhombohedral/data/";
if ~isfolder(default_root)
    default_root = "/Users/haoranyan/rg_master/rhombohedral/data/";
end
if ~isfolder(default_root)
    default_root = pwd;
end

iq_default = 4;
jq_default = 0;
FS = 14;

% =========================
% Choose x-axis
% =========================
x_mode = pick_xmode_();
fprintf("[x_mode] %s\n", x_mode);

% =========================
% Choose chi folder
% =========================
chi_root = uigetdir(default_root, ...
    "Select chi root folder containing chi*.txt recursively");

if isequal(chi_root,0)
    return;
end

chi_root = string(chi_root);

% =========================
% Choose iq, jq
% =========================
answer = inputdlg( ...
    {'iq:', 'jq:'}, ...
    'Choose q point', ...
    [1 35], ...
    {num2str(iq_default), num2str(jq_default)});

if isempty(answer)
    return;
end

iq_pick = str2double(answer{1});
jq_pick = str2double(answer{2});

if isnan(iq_pick) || isnan(jq_pick)
    error("Invalid iq/jq input.");
end

iq_pick = round(iq_pick);
jq_pick = round(jq_pick);

fprintf("[iq,jq] = (%d,%d)\n", iq_pick, jq_pick);

% =========================
% Load chi data
% =========================
[x, chi_re] = load_chi_vs_headerx_(chi_root, iq_pick, jq_pick, x_mode);

if isempty(x)
    error("No valid chi data found under: %s", chi_root);
end

[x, chi_re] = average_by_x_(x, chi_re);

% =========================
% Plot
% =========================
fig = figure('Color','w','Units','pixels','Position',[120 120 760 520]);
ax = axes(fig);
hold(ax,'on');
box(ax,'on');

% ❌ 不要 marker
plot(ax, x, chi_re, '-', ...
    'LineWidth', 2.0);

% 轴标签
if x_mode == "doping"
    xlabel(ax, 'doping (10^{12} cm^{-2})', ...
        'Interpreter','tex', 'FontSize', FS);
else
    xlabel(ax, '\mu (eV)', ...
        'Interpreter','tex', 'FontSize', FS);
end

ylabel(ax, 'Re(\chi)', ...
    'Interpreter','tex', 'FontSize', FS);

title(ax, sprintf('Re(\\chi), (iq,jq)=(%d,%d)', iq_pick, jq_pick), ...
    'Interpreter','tex', ...
    'FontSize', FS, ...
    'FontWeight','normal');

% ✅ 核心风格设置
set(ax, ...
    'FontSize', FS, ...
    'LineWidth', 1.2, ...
    'TickDir','in', ...     % ⭐ tick 朝里
    'Box','on');

% ❌ 去掉 grid
grid(ax,'off');

end

% ============================================================
% Helper: choose x-axis
% ============================================================
function x_mode = pick_xmode_()

choice = questdlg('Choose x-axis:', ...
    'x-axis', ...
    'doping', 'mu', 'doping');

if isempty(choice)
    x_mode = "doping";
else
    x_mode = string(choice);
end

end

% ============================================================
% Helper: load chi vs header x
% ============================================================
function [x, y] = load_chi_vs_headerx_(root_dir, iq_pick, jq_pick, x_mode)

root_dir = string(root_dir);

L = dir(fullfile(root_dir, "**", "chi*.txt"));

if isempty(L)
    x = [];
    y = [];
    return;
end

x = nan(numel(L),1);
y = nan(numel(L),1);

fprintf("[scan] found %d chi files\n", numel(L));

for k = 1:numel(L)

    fpath = string(fullfile(L(k).folder, L(k).name));

    H = parse_header_mu_doping_(fpath);

    if x_mode == "doping"
        xv = H.doping;
    else
        xv = H.mu;
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

    x(k) = xv;
    y(k) = M(id, C.Re);
end

m = isfinite(x) & isfinite(y);
x = x(m);
y = y(m);

fprintf("[kept] %d valid points\n", numel(x));

end

% ============================================================
% Helper: parse header mu / doping
% ============================================================
function meta = parse_header_mu_doping_(fpath)

meta = struct();
meta.doping = NaN;
meta.mu = NaN;

fid = fopen(fpath,'r');
if fid < 0
    return;
end

cleanupObj = onCleanup(@() fclose(fid));

num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';

for t = 1:5000
    ln = fgetl(fid);

    if ~ischar(ln)
        break;
    end

    s = strtrim(ln);

    if isempty(s)
        continue;
    end

    if s(1) ~= '#'
        break;
    end

    if ~isfinite(meta.doping)
        tok = regexp(s, ['^\#\s*doping\s*=\s*' num], ...
            'tokens', 'once');
        if ~isempty(tok)
            meta.doping = str2double(tok{1});
        end
    end

    if ~isfinite(meta.mu)
        tok = regexp(s, ['^\#\s*mu\s*=\s*' num], ...
            'tokens', 'once');

        if isempty(tok)
            tok = regexp(s, ['^\#\s*EF\s*=\s*' num], ...
                'tokens', 'once');
        end

        if ~isempty(tok)
            meta.mu = str2double(tok{1});
        end
    end

    if isfinite(meta.doping) && isfinite(meta.mu)
        break;
    end
end

end

% ============================================================
% Helper: read numeric data, skipping # lines
% ============================================================
function M = read_numeric_skiphash_(fpath)

fid = fopen(fpath,'r');

if fid < 0
    M = [];
    return;
end

cleanupObj = onCleanup(@() fclose(fid));

rows = {};

while true
    ln = fgetl(fid);

    if ~ischar(ln)
        break;
    end

    if isempty(strtrim(ln))
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

lastFinite = find(any(isfinite(M),1), 1, 'last');

if ~isempty(lastFinite)
    M = M(:,1:lastFinite);
end

end

% ============================================================
% Helper: detect chi columns
% ============================================================
function C = detect_cols_chi_(M)

C = struct();

% Format A:
% iq jq qx qy chi_real chi_imag nKpair nK
% Example: 12 -6 0 -0.0153 0.0730 ...
if size(M,2) >= 6
    C.iq = 1;
    C.jq = 2;
    C.Re = 5;
else
    error("Chi numeric data has too few columns.");
end

end

% ============================================================
% Helper: average duplicate x
% ============================================================
function [x_u, y_u] = average_by_x_(x, y)

if isempty(x)
    x_u = x;
    y_u = y;
    return;
end

xq = round(x, 12);

[xs, ord] = sort(xq);
ys = y(ord);

[ux, ~, ic] = unique(xs);

y_u = accumarray(ic, ys, [], @(v) mean(v,'omitnan'));
x_u = ux;

end