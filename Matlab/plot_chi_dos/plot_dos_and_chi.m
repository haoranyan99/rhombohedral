function plot_dos_and_chi()
% Plot DOS and Re(chi) on the same figure vs x-axis = doping or mu.
% Key rule: x for chi MUST come from each file HEADER (not folder name).
% DOS: prefers x from DOS numeric columns; if missing, tries header scalar.

default_root = "/Users/haoranyan/rg_master/data/";

x_mode = pick_xmode_();  % "doping" or "mu"
fprintf("[x_mode] %s\n", x_mode);

% -----------------------------
% default iq/jq
% -----------------------------
iq0 = -267;
jq0 = -267;
% iq0 = 1;
% jq0 = 1;

% =========================
% 1) DOS file (single)
% =========================
start_dos = default_root;
if ~isfolder(start_dos)
    start_dos = string(pwd);
end

[fname, fpath] = uigetfile({'*.txt','DOS data (*.txt)'}, 'Select DOS file', start_dos);
if isequal(fname,0), return; end
dos_path = fullfile(fpath, fname);
[~, dos_base, ~] = fileparts(dos_path);

raw = readmatrix(dos_path, "FileType","text", "CommentStyle","#");
if isempty(raw) || size(raw,2) < 3
    error("Bad DOS file: need numeric block with >=3 columns.");
end

cols_line = parse_columns_line_(dos_path); % "" if not found
meta_dos  = parse_header_mu_doping_(dos_path); % may be NaN for missing

[x_dos, DOS] = extract_x_and_dos_(raw, cols_line, x_mode, meta_dos);

m = isfinite(x_dos) & isfinite(DOS);
x_dos = x_dos(m);
DOS   = DOS(m);

[x_dos, DOS] = average_by_x_(x_dos, DOS);

% =========================
% 2) chi folder (many files)
% =========================
start_chi = default_root;
if ~isfolder(start_chi)
    start_chi = string(fpath);
end

chi_root = uigetdir(start_chi, "Select chi root folder containing chi*.txt (recursive)");
if isequal(chi_root,0), return; end
chi_root = string(chi_root);

iq_pick = iq0;
jq_pick = jq0;

fprintf("[iq,jq] = (%d, %d)\n", iq_pick, jq_pick);

[x_chi, chi_re] = load_chi_vs_headerx_(chi_root, iq_pick, jq_pick, x_mode);
if isempty(x_chi)
    error("No valid chi points found under: %s (x_mode=%s). Check headers contain that x.", chi_root, x_mode);
end
[x_chi, chi_re] = average_by_x_(x_chi, chi_re);

% =========================
% 3) overlap range
% =========================
xmin = max(min(x_dos), min(x_chi));
xmax = min(max(x_dos), max(x_chi));
if ~(xmax > xmin)
    error("No overlap in %s range between DOS and chi.", x_mode);
end

md = (x_dos >= xmin) & (x_dos <= xmax);
mc = (x_chi >= xmin) & (x_chi <= xmax);

xx_d = x_dos(md); yy_d = DOS(md);
xx_c = x_chi(mc); yy_c = chi_re(mc);

% =========================
% 4) plot
% =========================
FS = 14;
fig = figure('Color','w','Units','pixels','Position',[120 120 860 540]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'FontSize',FS,'TickDir','out','LineWidth',1.0);

yyaxis(ax,'left');
plot(ax, xx_d, yy_d, '-','LineWidth',2);
ylabel(ax,'DOS (arb.)','Interpreter','none','FontSize',FS);

yyaxis(ax,'right');
plot(ax, xx_c, yy_c, 'o-','LineWidth',1.8,'MarkerSize',5);
ylabel(ax,'Re(\chi) (arb.)','Interpreter','none','FontSize',FS);

if x_mode=="doping"
    xlabel(ax,'doping (1e12 cm^{-2})','Interpreter','none','FontSize',FS);
else
    xlabel(ax,'\mu (eV)','Interpreter','tex','FontSize',FS);
end

title(ax, sprintf('%s | (iq,jq)=(%d,%d) | overlap=[%.6g, %.6g]', ...
    dos_base, iq_pick, jq_pick, xmin, xmax), 'Interpreter','none','FontWeight','normal');

legend(ax, {'DOS','Re(\chi)'}, 'Location','best');
end

% ============================================================
% ---- helpers (minimal, standalone) ----
% ============================================================

function x_mode = pick_xmode_()
c = questdlg('Choose x-axis:', 'x-axis', 'doping','mu','doping');
if isempty(c), x_mode = "doping"; else, x_mode = string(c); end
end

function cols_line = parse_columns_line_(in_path)
cols_line = "";
fid = fopen(in_path,'r');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    s = strtrim(ln);
    if ~startsWith(s,"#"), break; end
    s2 = strtrim(erase(s,"#"));
    if startsWith(lower(s2), "columns:")
        cols_line = string(s2);
        return;
    end
end
end

function meta = parse_header_mu_doping_(fpath)
% Look for:
%   # doping = ...
%   # mu = ...  or # EF = ...
meta = struct('doping',NaN,'mu',NaN);

fid = fopen(fpath,'r');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';

for t = 1:4000
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    s = strtrim(ln);
    if isempty(s), continue; end
    if s(1) ~= '#', break; end

    if ~isfinite(meta.doping)
        tok = regexp(s, ['^\#\s*doping\s*=\s*' num], 'tokens','once');
        if ~isempty(tok), meta.doping = str2double(tok{1}); end
    end

    if ~isfinite(meta.mu)
        tok = regexp(s, ['^\#\s*mu\s*=\s*' num], 'tokens','once');
        if isempty(tok)
            tok = regexp(s, ['^\#\s*EF\s*=\s*' num], 'tokens','once');
        end
        if ~isempty(tok), meta.mu = str2double(tok{1}); end
    end

    if isfinite(meta.doping) && isfinite(meta.mu), break; end
end
end

function [x, DOS] = extract_x_and_dos_(raw, cols_line, x_mode, meta_dos)
% Prefer x from numeric columns. If missing, try header scalar (meta_dos).
ncol = size(raw,2);
DOS = raw(:,end);
x = nan(size(DOS));

if strlength(cols_line) > 0
    s = char(cols_line);
    s = strrep(s, "columns:", "");
    s = strtrim(s);
    s = regexprep(s, "\s+", " ");
    toks = strsplit(s, " ");

    % DOS index
    idx_dos = find(strcmpi(toks,"DOS") | strcmpi(toks,"dos"), 1);
    if isempty(idx_dos), idx_dos = ncol; end
    DOS = raw(:, min(idx_dos,ncol));

    if x_mode=="doping"
        idx_x = find(contains(lower(toks),"doping"), 1);
        if ~isempty(idx_x) && idx_x <= ncol
            x = raw(:, idx_x);
            return;
        end
        % fallback: header scalar doping (constant across rows)
        if isfinite(meta_dos.doping)
            x = meta_dos.doping * ones(size(DOS));
            return;
        end
        error("DOS: cannot find doping column in numeric table nor 'doping=' in header.");
    else
        % mu
        idx_x = find(contains(lower(toks),"e(ev)") | contains(lower(toks),"ef"), 1);
        if ~isempty(idx_x) && idx_x <= ncol
            x = raw(:, idx_x);
            return;
        end
        if isfinite(meta_dos.mu)
            x = meta_dos.mu * ones(size(DOS));
            return;
        end
        error("DOS: cannot find mu/EF column in numeric table nor 'mu='/'EF=' in header.");
    end
end

% no columns line: try simple fallback
if x_mode=="doping"
    % common format: i E filling doping DOS  => doping is col4
    if ncol >= 4
        x = raw(:,4);
        return;
    end
    if isfinite(meta_dos.doping)
        x = meta_dos.doping * ones(size(DOS));
        return;
    end
    error("DOS: no '# columns:' and cannot infer doping column; also no header doping.");
else
    if isfinite(meta_dos.mu)
        x = meta_dos.mu * ones(size(DOS));
        return;
    end
    error("DOS: no '# columns:' and cannot infer mu column; also no header mu/EF.");
end
end

function [x_u, y_u] = average_by_x_(x, y)
if isempty(x), x_u=x; y_u=y; return; end
xq = round(x, 12);
[xs, ord] = sort(xq);
ys = y(ord);
[ux, ~, ic] = unique(xs);
y_u = accumarray(ic, ys, [], @(v) mean(v,'omitnan'));
x_u = ux;
end

function [x, y] = load_chi_vs_headerx_(root_dir, iq_pick, jq_pick, x_mode)
% Scan chi*.txt under root_dir recursively.
% For each file:
%   - read header: T, doping, mu/EF (we only need x_mode one)
%   - read numeric: find (iq,jq) and take Re(chi)
% Return vectors x (from header) and y=Re(chi).

root_dir = string(root_dir);
L = dir(fullfile(root_dir, "**", "chi*.txt"));
if isempty(L), x=[]; y=[]; return; end

x = nan(numel(L),1);
y = nan(numel(L),1);

for k = 1:numel(L)
    fpath = string(fullfile(L(k).folder, L(k).name));

    H = parse_header_mu_doping_(fpath);
    if x_mode=="doping"
        xv = H.doping;
    else
        xv = H.mu;
    end
    if ~isfinite(xv), continue; end

    M = read_numeric_skiphash_(fpath);
    if isempty(M) || size(M,2) < 6, continue; end

    Cc = detect_cols_chi_(M);
    id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
    if isempty(id), continue; end

    x(k) = xv;
    y(k) = M(id, Cc.Re);
end

m = isfinite(x) & isfinite(y);
x = x(m); y = y(m);
end

function M = read_numeric_skiphash_(fpath)
fid = fopen(fpath,'r');
if fid < 0, M=[]; return; end
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
if isempty(rows), M=[]; return; end
ncol = max(cellfun(@numel, rows));
M = nan(numel(rows), ncol);
for i = 1:numel(rows)
    v = rows{i};
    M(i,1:numel(v)) = v;
end
lastFinite = find(any(isfinite(M),1), 1, 'last');
if ~isempty(lastFinite), M = M(:,1:lastFinite); end
end

function C = detect_cols_chi_(M)
% 7: iq jq qx qy Re Im ...
% 8+: idx iq jq qx qy Re Im ...
ncol = size(M,2);
C = struct('iq',1,'jq',2,'Re',5);
if ncol >= 8
    C.iq = 2; C.jq = 3; C.Re = 6;
end
end