function plot_dos_and_chi_multi_pick()
% Plot DOS and multiple Re(chi) curves on the same figure vs x-axis = doping or mu.
% Workflow:
%   1) choose DOS file once
%   2) repeatedly:
%        - input (iq,jq)
%        - choose chi folder for this (iq,jq)
%        - load and append this chi curve
%   3) all chi curves share the same right axis
%
% x for chi MUST come from each file HEADER (not folder name).

default_root = "/Users/haoranyan/rg_master/data/";

x_mode = pick_xmode_();  % "doping" or "mu"
fprintf("[x_mode] %s\n", x_mode);

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

cols_line = parse_columns_line_(dos_path);
meta_dos  = parse_header_mu_doping_(dos_path);

[x_dos, DOS] = extract_x_and_dos_(raw, cols_line, x_mode, meta_dos);

m = isfinite(x_dos) & isfinite(DOS);
x_dos = x_dos(m);
DOS   = DOS(m);

[x_dos, DOS] = average_by_x_(x_dos, DOS);

% =========================
% 2) repeatedly choose chi folders for different (iq,jq)
% =========================
chi_sets = struct('iq',{},'jq',{},'root',{},'x',{},'y',{});

iq0 = -131;
jq0 = -131;

iq_pick = iq0;
jq_pick = jq0;

while true
    prompt = sprintf(['Current (iq,jq) = (%d, %d)\n\n' ...
        'Select chi root folder for this (iq,jq).\n' ...
        'Cancel = stop adding chi curves.'], iq_pick, jq_pick);

    chi_root = uigetdir(default_root, prompt);
    if isequal(chi_root,0)
        break;
    end
    chi_root = string(chi_root);

    fprintf("[iq,jq] = (%d, %d), folder = %s\n", iq_pick, jq_pick, chi_root);

    [x_chi, chi_re] = load_chi_vs_headerx_(chi_root, iq_pick, jq_pick, x_mode);
    if isempty(x_chi)
        warning("No valid chi points found under: %s for (iq,jq)=(%d,%d)", ...
            chi_root, iq_pick, jq_pick);
    else
        [x_chi, chi_re] = average_by_x_(x_chi, chi_re);

        chi_sets(end+1).iq   = iq_pick; %#ok<AGROW>
        chi_sets(end).jq     = jq_pick;
        chi_sets(end).root   = chi_root;
        chi_sets(end).x      = x_chi;
        chi_sets(end).y      = chi_re;
    end

    ans_ij = inputdlg( ...
        {'Next iq:','Next jq:'}, ...
        'Input next (iq,jq)', ...
        [1 35], ...
        {num2str(iq_pick), num2str(jq_pick)} );

    if isempty(ans_ij)
        break;
    end

    iq_new = str2double(ans_ij{1});
    jq_new = str2double(ans_ij{2});

    if ~isfinite(iq_new) || ~isfinite(jq_new)
        warning('Invalid next (iq,jq). Stop adding chi curves.');
        break;
    end

    iq_pick = round(iq_new);
    jq_pick = round(jq_new);
end

if isempty(chi_sets)
    error("No valid chi curves were added.");
end

% =========================
% 3) overlap range
% =========================
xmin = min(x_dos);
xmax = max(x_dos);

for k = 1:numel(chi_sets)
    xmin = max(xmin, min(chi_sets(k).x));
    xmax = min(xmax, max(chi_sets(k).x));
end

if ~(xmax > xmin)
    error("No overlap in %s range between DOS and chi curves.", x_mode);
end

md = (x_dos >= xmin) & (x_dos <= xmax);
xx_d = x_dos(md);
yy_d = DOS(md);

for k = 1:numel(chi_sets)
    mc = (chi_sets(k).x >= xmin) & (chi_sets(k).x <= xmax);
    chi_sets(k).x = chi_sets(k).x(mc);
    chi_sets(k).y = chi_sets(k).y(mc);
end

% =========================
% 4) plot
% =========================
FS = 14;
fig = figure('Color','w','Units','pixels','Position',[120 120 900 560]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'FontSize',FS,'TickDir','out','LineWidth',1.0);

yyaxis(ax,'left');
h_dos = plot(ax, xx_d, yy_d, '-k', 'LineWidth', 2);
ylabel(ax,'DOS (arb.)','Interpreter','none','FontSize',FS);

yyaxis(ax,'right');
cc = lines(numel(chi_sets));
h_chi = gobjects(numel(chi_sets),1);
leg = cell(numel(chi_sets)+1,1);
leg{1} = 'DOS';

for k = 1:numel(chi_sets)
    h_chi(k) = plot(ax, chi_sets(k).x, chi_sets(k).y, 'o-', ...
        'LineWidth',1.8, 'MarkerSize',5, 'Color',cc(k,:));
    [~, folder_name] = fileparts(char(chi_sets(k).root));
    leg{k+1} = sprintf('Re(chi) (%d,%d) | %s', ...
        chi_sets(k).iq, chi_sets(k).jq, folder_name);
end
ylabel(ax,'Re(\chi) (arb.)','Interpreter','none','FontSize',FS);

if x_mode=="doping"
    xlabel(ax,'doping (1e12 cm^{-2})','Interpreter','none','FontSize',FS);
else
    xlabel(ax,'\mu (eV)','Interpreter','tex','FontSize',FS);
end

title(ax, sprintf('%s | overlap=[%.6g, %.6g]', ...
    dos_base, xmin, xmax), ...
    'Interpreter','none','FontWeight','normal');

legend(ax, [h_dos; h_chi], leg, 'Location','best', 'Interpreter','none');
end

% ============================================================
% ---- helpers ----
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
ncol = size(raw,2);
DOS = raw(:,end);
x = nan(size(DOS));

if strlength(cols_line) > 0
    s = char(cols_line);
    s = strrep(s, "columns:", "");
    s = strtrim(s);
    s = regexprep(s, "\s+", " ");
    toks = strsplit(s, " ");

    idx_dos = find(strcmpi(toks,"DOS") | strcmpi(toks,"dos"), 1);
    if isempty(idx_dos), idx_dos = ncol; end
    DOS = raw(:, min(idx_dos,ncol));

    if x_mode=="doping"
        idx_x = find(contains(lower(toks),"doping"), 1);
        if ~isempty(idx_x) && idx_x <= ncol
            x = raw(:, idx_x);
            return;
        end
        if isfinite(meta_dos.doping)
            x = meta_dos.doping * ones(size(DOS));
            return;
        end
        error("DOS: cannot find doping column in numeric table nor 'doping=' in header.");
    else
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

if x_mode=="doping"
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
ncol = size(M,2);
C = struct('iq',1,'jq',2,'Re',5);
if ncol >= 8
    C.iq = 2; C.jq = 3; C.Re = 6;
end
end