function out = plot_chi_max_vs_polar_mu()
% plot_chi_max_vs_polar_mu
% ------------------------------------------------------------
% For each polar_mu (meV): collect chi(iq,jq) at ONE T (T_target±T_tol)
% under filters (mu_range, dop_range), then compute:
%   chi_max(polar) = max(chi)  OR  max(abs(chi))
% Plot chi_max vs polar_mu.
% ------------------------------------------------------------

out = struct();
out.png = "";

% =========================
% USER SETTINGS (EDIT HERE)
% =========================
D_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_200\D0.067";

T_target = 4;        % K
T_tol    = 1e-6;     % K

iq_pick = -133;
jq_pick = -133;

polar_range = [];            % meV, [] => all
mu_range    = [];            % eV,  [] => all
dop_range   = [];            % [] => all

max_mode = "max";            % "max" or "maxabs"

save_png = true;
out_dir  = fullfile(string(D_root), "plot_chiMax_vs_polar");
if save_png && ~exist(out_dir,"dir"), mkdir(out_dir); end
% =========================

assert(isfolder(D_root), "D_root not found: %s", D_root);

fprintf("Root:\n  %s\n", D_root);
fprintf("Pick (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);
fprintf("T_target=%.12g ± %.3g\n", T_target, T_tol);
fprintf("filters: polar=%s | mu=%s | dop=%s\n", range_str(polar_range), range_str(mu_range), range_str(dop_range));
fprintf("max_mode = %s\n", max_mode);

T_lo = T_target - T_tol;
T_hi = T_target + T_tol;

% -----------------------------
% scan chi files
% -----------------------------
files = dir(fullfile(D_root, "**", "chi*.txt"));
files = files(~[files.isdir]);
Nraw = numel(files);
if Nraw==0, error("No chi*.txt found under: %s", D_root); end
fprintf("[scan] %d chi files\n", Nraw);

% -----------------------------
% collect by polar_mu
%   rec(key).polar0
%   rec(key).chi_all
% -----------------------------
rec = containers.Map('KeyType','char','ValueType','any');

n_pass=0; n_skip_noT=0; n_skip_filter=0; n_skip_read=0; n_skip_nomatch=0;

for k = 1:Nraw
    p = string(fullfile(files(k).folder, files(k).name));

    m = parse_header_T_doping_mu(p);
    if ~isfinite(m.T)
        n_skip_noT = n_skip_noT + 1;
        continue;
    end

    if ~(m.T >= T_lo && m.T <= T_hi)
        n_skip_filter = n_skip_filter + 1;
        continue;
    end

    polar = parse_polar_from_path(p);

    if ~pass_range_AND(polar, polar_range) || ...
       ~pass_range_AND(m.mu, mu_range) || ...
       ~pass_range_AND(m.doping, dop_range)
        n_skip_filter = n_skip_filter + 1;
        continue;
    end

    [iq_col, jq_col, re_col] = read_chi_table_cols(p);
    if isempty(iq_col)
        n_skip_read = n_skip_read + 1;
        continue;
    end

    tf = (iq_col(:) == iq_pick) & (jq_col(:) == jq_pick);
    if ~any(tf)
        n_skip_nomatch = n_skip_nomatch + 1;
        continue;
    end

    chi = double(re_col(find(tf,1,'first')));

    key = make_polar_key(polar);
    add_chi_by_polar(rec, key, polar, chi);

    n_pass = n_pass + 1;
end

fprintf("[done] passed=%d | skip(noT)=%d | skip(filter)=%d | skip(read)=%d | skip(no_ij)=%d | polars=%d\n", ...
    n_pass, n_skip_noT, n_skip_filter, n_skip_read, n_skip_nomatch, rec.Count);

if rec.Count==0
    error("No data left after filters. Relax T_tol / ranges / iq,jq.");
end

% -----------------------------
% compute chi_max per polar
% -----------------------------
keys = rec.keys;
nP = numel(keys);
polar_list = nan(nP,1);
chi_max    = nan(nP,1);
npts       = nan(nP,1);

for i = 1:nP
    S = rec(keys{i});
    polar_list(i) = S.polar0;
    v = S.chi_all(:);
    v = v(isfinite(v));
    npts(i) = numel(v);

    if isempty(v)
        chi_max(i) = NaN;
    else
        switch string(max_mode)
            case "max"
                chi_max(i) = max(v);
            case "maxabs"
                chi_max(i) = max(abs(v));
            otherwise
                error('max_mode must be "max" or "maxabs".');
        end
    end
end

% sort by polar
[polar_list, ord] = sort(polar_list);
chi_max = chi_max(ord);
npts = npts(ord);

ok = isfinite(polar_list) & isfinite(chi_max);
if nnz(ok)<1
    error("No finite polar/chi_max to plot.");
end

% -----------------------------
% plot
% -----------------------------
FS = 14;
fig = figure("Color","w","Units","pixels","Position",[180 120 820 560], ...
    "Name","chi_max vs polar_mu");
ax = axes(fig); hold(ax,"on"); box(ax,"on"); grid(ax,"on");

plot(ax, polar_list(ok), chi_max(ok), "o-","LineWidth",2,"MarkerSize",6);

xlabel(ax, "$\mathrm{polar\_mu}\ (\mathrm{meV})$", "Interpreter","latex");
switch string(max_mode)
    case "max"
        ylabel(ax, sprintf("$\\max\\,\\mathrm{Re}\\,\\chi(%d,%d)$", iq_pick, jq_pick), "Interpreter","latex");
        modeStr = "max(Re chi)";
    case "maxabs"
        ylabel(ax, sprintf("$\\max\\,|\\mathrm{Re}\\,\\chi(%d,%d)|$", iq_pick, jq_pick), "Interpreter","latex");
        modeStr = "max(|Re chi|)";
end

title(ax, sprintf("$T=%.6g\\,\\mathrm{K}$, $q=(%d,%d)$, %s", T_target, iq_pick, jq_pick, modeStr), ...
    "Interpreter","latex","FontWeight","normal");

set(ax,"FontSize",FS,"TickLabelInterpreter","latex","LineWidth",1.0,"TickDir","out");

% optional: annotate point counts
for i = 1:numel(polar_list)
    if ok(i)
        text(ax, polar_list(i), chi_max(i), sprintf("  n=%d", npts(i)), ...
            "FontSize",FS-4,"Interpreter","none");
    end
end

% save
if save_png
    out.png = fullfile(out_dir, sprintf("chiMax_vs_polar_T%.6g_i%d_j%d_%s.png", ...
        T_target, iq_pick, jq_pick, max_mode));
    exportgraphics(fig, out.png, "Resolution", 250);
    fprintf("[saved] %s\n", out.png);
end

out.polar_meV = polar_list;
out.chi_max = chi_max;
out.n_points = npts;
out.max_mode = max_mode;
out.T_target = T_target;
out.iq_pick = iq_pick;
out.jq_pick = jq_pick;
out.filters = struct("polar_range",polar_range,"mu_range",mu_range,"dop_range",dop_range);

end

% ============================================================
% helpers (same style as yours)
% ============================================================
function s = range_str(r)
    if isempty(r)
        s = "[]";
    else
        s = sprintf("[%.12g, %.12g]", min(r), max(r));
    end
end

function tf = pass_range_AND(v, range)
    if isempty(range), tf = true; return; end
    if ~isfinite(v), tf = false; return; end
    lo = min(range); hi = max(range);
    tf = (v >= lo) && (v <= hi);
end

function polar = parse_polar_from_path(p)
    polar = NaN;
    tok = regexp(p, "[\\/](?:polar_meV|polar)\s*([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)", "tokens", "once");
    if ~isempty(tok)
        polar = str2double(tok{1});
    end
end

function m = parse_header_T_doping_mu(in_path)
    m = struct("T",NaN,"doping",NaN,"mu",NaN);

    fid = fopen(in_path,'r');
    if fid < 0, return; end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end
        s0 = strtrim(tline);
        if ~startsWith(s0, "#"), break; end
        s = strtrim(erase(s0, "#"));

        if ~isfinite(m.T)
            tok = regexp(s, "(?:^|\s)T(?:_K)?\s*=\s*" + num, "tokens", "once");
            if ~isempty(tok), m.T = str2double(tok{1}); end
        end
        if ~isfinite(m.doping)
            tok = regexp(s, "doping\s*=\s*" + num, "tokens", "once");
            if ~isempty(tok), m.doping = str2double(tok{1}); end
        end
        if ~isfinite(m.mu)
            tok = regexp(s, "(?:^|\s)(?:mu|mu_eV|EF|EF_eV)\s*=\s*" + num, "tokens", "once");
            if ~isempty(tok), m.mu = str2double(tok{1}); end
        end
    end
end

function [iq, jq, re] = read_chi_table_cols(in_path)
    iq = []; jq = []; re = [];
    try
        raw = readmatrix(in_path, "FileType","text", "CommentStyle","#");
    catch
        return;
    end
    if isempty(raw) || size(raw,2) < 6
        return;
    end

    if size(raw,2) >= 8
        iq = raw(:,2);
        jq = raw(:,3);
        re = raw(:,6);
    else
        iq = raw(:,1);
        jq = raw(:,2);
        re = raw(:,5);
    end
end

function key = make_polar_key(polar)
    if ~isfinite(polar)
        key = "polar=nan";
    else
        key = sprintf("polar=%.12g", polar);
    end
end

function add_chi_by_polar(M, key, polar0, chi)
    if ~isKey(M, key)
        S = struct();
        S.polar0 = polar0;
        S.chi_all = chi;
        M(key) = S;
    else
        S = M(key);
        S.chi_all(end+1,1) = chi;
        M(key) = S;
    end
end