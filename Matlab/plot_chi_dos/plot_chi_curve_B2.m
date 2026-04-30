function out = plot_chi_curve_B2()
% plot_chi_Bcurves_oneT_onePoint_noUI
% Fixed: one temperature (T_target) + one (iq,jq) point.
% Curves: different polar_mu as different curves.
%
% Output:
%   - FigMu:  Re chi(iq,jq) vs mu      (curves grouped by polar_mu)
%   - FigDop: Re chi(iq,jq) vs doping  (curves grouped by polar_mu)

out = struct();  
out.png_mu  = "";
out.png_dop = "";

% ============================================================
% USER SETTINGS (EDIT HERE)
% ============================================================

% D_root = "/Users/haoranyan/Library/CloudStorage/OneDrive-Emory/Rhombohedral_SC/rhombohedral_project/data/chi_sk_mu_200/D0.084";
default_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_200\D0.084";
if ~isfolder(default_root)
    warning("Default folder not found: %s\nFallback to pwd.", default_root);
    default_root = string(pwd);
end

D_root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
if isequal(D_root, 0), error('User cancelled.'); end
D_root = string(D_root);
fprintf("Root folder:\n  %s\n", D_root);

% ---- pick ONE temperature ----
T_target = 6.5;        % K
T_tol    = 1e-6;     % K tolerance

% ---- pick ONE (iq,jq) ----
iq_pick = 3;
jq_pick = -3;

% ---- AND filters; [] => disable ----
polar_range = 0:1:6;       % meV; [] => all polar
mu_range    = [];    % eV;  [] => all mu
dop_range   = [-2.5:0];           % [] => all doping

% ---- IMPORTANT: speed switch ----
% If true: ONLY scan chi*.txt under folders named polar_meV*
% This avoids the painfully slow dir("**") across the whole OneDrive tree.
only_polar_meV_folders = true;

% ---- output ----
save_png = true;
out_dir  = fullfile(string(D_root), "plot_Bcurves_oneT_onePoint");
if save_png && ~exist(out_dir,"dir"), mkdir(out_dir); end

% ---- style ----
FS = 14;
use_colorbar = false; % 你想要colorbar就改成 true

% ============================================================

assert(isfolder(D_root), "D_root not found: %s", D_root);

fprintf("Root D folder:\n  %s\n", D_root);
fprintf("Pick (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);
fprintf("Pick T_target=%.12g K, tol=%.3g K => range=[%.12g, %.12g]\n", ...
    T_target, T_tol, T_target-T_tol, T_target+T_tol);
fprintf("[filter AND]\n");
fprintf("  polar_range = %s\n", range_str(polar_range));
fprintf("  mu_range    = %s\n", range_str(mu_range));
fprintf("  dop_range   = %s\n", range_str(dop_range));

% -----------------------------
% scan chi files (FAST MODE)
% -----------------------------
t_scan = tic;
files = [];

if only_polar_meV_folders
    polarDirs = dir(fullfile(D_root, "**", "polar_meV*"));
    polarDirs = polarDirs([polarDirs.isdir]);

    fprintf("[scan] polar_meV* folders found: %d\n", numel(polarDirs));
    if isempty(polarDirs)
        error("No polar_meV* folders found under: %s", D_root);
    end

    % Collect chi*.txt under these polar folders (still recursive, but much smaller subtrees)
    for i = 1:numel(polarDirs)
        base = fullfile(polarDirs(i).folder, polarDirs(i).name);
        ff = dir(fullfile(base, "**", "chi*.txt"));
        ff = ff(~[ff.isdir]);
        if ~isempty(ff)
            files = [files; ff]; %#ok<AGROW>
        end
        if mod(i, 20) == 0 || i == numel(polarDirs)
            fprintf("  scanning polar folders: %d/%d | chi files collected so far: %d | elapsed %.1fs\n", ...
                i, numel(polarDirs), numel(files), toc(t_scan));
            drawnow;
        end
    end
else
    % Original behavior (can be very slow on OneDrive trees)
    files = dir(fullfile(D_root, "**", "chi*.txt"));
    files = files(~[files.isdir]);
end

Nraw = numel(files);
if Nraw == 0
    error("No chi*.txt found under: %s (mode only_polar_meV_folders=%d)", D_root, only_polar_meV_folders);
end
fprintf("[scan] found %d chi files | scan time = %.1fs\n", Nraw, toc(t_scan));

% -----------------------------
% collect by polar_mu:
%   rec(polar_key) stores arrays: mu, dop, chi, paths, polar0
% -----------------------------
rec = containers.Map('KeyType','char','ValueType','any');

t0 = tic;
t_last = 0;
print_every = 0.5;

n_pass = 0;
n_skip_noT = 0;
n_skip_filter = 0;
n_skip_read = 0;
n_skip_nomatch = 0;

T_lo = T_target - T_tol;
T_hi = T_target + T_tol;

for k = 1:Nraw
    p = string(fullfile(files(k).folder, files(k).name));

    m = parse_header_T_doping_mu(p);
    if ~isfinite(m.T)
        n_skip_noT = n_skip_noT + 1;
        progress_print();
        continue;
    end

    % fixed T filter
    if ~(m.T >= T_lo && m.T <= T_hi)
        n_skip_filter = n_skip_filter + 1;
        progress_print();
        continue;
    end

    polar = parse_polar_from_path(p);

    % strict AND filters
    if ~pass_range_AND(polar, polar_range) || ...
       ~pass_range_AND(m.mu, mu_range) || ...
       ~pass_range_AND(m.doping, dop_range)
        n_skip_filter = n_skip_filter + 1;
        progress_print();
        continue;
    end

    % read numeric table once
    [iq_col, jq_col, re_col] = read_chi_table_cols(p);
    if isempty(iq_col)
        n_skip_read = n_skip_read + 1;
        progress_print();
        continue;
    end

    % find this (iq_pick,jq_pick)
    tf = (iq_col(:) == iq_pick) & (jq_col(:) == jq_pick);
    if ~any(tf)
        n_skip_nomatch = n_skip_nomatch + 1;
        progress_print();
        continue;
    end
    chi = double(re_col(find(tf,1,'first')));

    key = make_polar_key(polar);
    add_point_by_polar(rec, key, polar, m.mu, m.doping, chi, p);
    n_pass = n_pass + 1;

    progress_print();
end

fprintf("[done] passed=%d | skipped(noT)=%d | skipped(filter)=%d | skipped(read)=%d | skipped(no_ij)=%d | curves(polar)=%d | time=%.1fs\n", ...
    n_pass, n_skip_noT, n_skip_filter, n_skip_read, n_skip_nomatch, rec.Count, toc(t0));

if rec.Count == 0
    error("No curves produced: check T_target/T_tol, filters, or (iq,jq).");
end

% -----------------------------
% color mapping by polar_mu (global)
% -----------------------------
keys = rec.keys;
polar_list = nan(numel(keys),1);
for i = 1:numel(keys)
    S = rec(keys{i});
    polar_list(i) = S.polar0;
end

Pmin = min(polar_list(isfinite(polar_list)));
Pmax = max(polar_list(isfinite(polar_list)));
if ~isfinite(Pmin) || ~isfinite(Pmax)
    error("No finite polar_mu found in curves.");
end
if Pmin == Pmax, Pmax = Pmin + 1e-12; end

% sort curves by polar for nicer legend order
[polar_sorted, ord] = sort(polar_list); %#ok<ASGLU>
keys = keys(ord);

cmap = turbo(256);

set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
set(groot, "defaultTextInterpreter", "latex");

% ============================================================
% FigDop: chi vs doping
% ============================================================
figD = figure("Color","w","Units","pixels","Position",[120 120 1200 750], ...
    "Name","chi(iq,jq) vs doping (polar curves)");
axD = axes(figD,"Position",[0.08 0.12 0.82 0.82]); hold(axD,"on");
colormap(axD, cmap);

nplotted_d = 0;
for i = 1:numel(keys)
    S = rec(keys{i});
    [x, y] = collapse_same_x(S.dop, S.chi);
    if numel(x) < 2, continue; end

    pval = S.polar0;
    ci = round((pval - Pmin)/(Pmax - Pmin) * 255) + 1;
    ci = max(min(ci,256),1);

    plot(axD, x, y, "-", ...
        "LineWidth", 1.8, ...
        "Color", cmap(ci,:), ...
        "DisplayName", sprintf("polar = %.6g\\,meV", pval));
    nplotted_d = nplotted_d + 1;
end

xlabel(axD,"doping (header units)","FontSize",FS);
ylabel(axD,sprintf("Re \\chi(%d,%d)", iq_pick, jq_pick),"FontSize",FS);
title(axD, sprintf("Re \\chi vs doping | T=%.6g K (tol %.1g) | curves=%d", ...
    T_target, T_tol, nplotted_d), "FontSize",FS, "FontWeight","normal");
grid(axD,"on"); box(axD,"on");
set(axD,"FontSize",FS,"LineWidth",1.0,"TickDir","out");

if use_colorbar
    caxis(axD,[Pmin Pmax]);
    cb = colorbar(axD);
    cb.Label.String = "polar\\_mu (meV)";
    cb.Label.Interpreter = "latex";
end

lg = legend(axD, "Location", "eastoutside");
lg.Title.String = "polar\\_mu";
lg.Title.Interpreter = "latex";
lg.FontSize = FS-2;

if save_png
    out.png_dop = fullfile(out_dir, sprintf("chi_vs_doping_T%.6g_i%d_j%d.png", T_target, iq_pick, jq_pick));
    exportgraphics(figD, out.png_dop, "Resolution", 220);
    fprintf("[saved] %s\n", out.png_dop);
end

% ============================================================
% FigMu: chi vs mu
% ============================================================
figM = figure("Color","w","Units","pixels","Position",[160 160 1200 750], ...
    "Name","chi(iq,jq) vs mu (polar curves)");
axM = axes(figM,"Position",[0.08 0.12 0.82 0.82]); hold(axM,"on");
colormap(axM, cmap);

nplotted_m = 0;
for i = 1:numel(keys)
    S = rec(keys{i});
    [x, y] = collapse_same_x(S.mu, S.chi);
    if numel(x) < 2, continue; end

    pval = S.polar0;
    ci = round((pval - Pmin)/(Pmax - Pmin) * 255) + 1;
    ci = max(min(ci,256),1);

    plot(axM, x, y, "-", ...
        "LineWidth", 1.8, ...
        "Color", cmap(ci,:), ...
        "DisplayName", sprintf("polar = %.6g\\,meV", pval));
    nplotted_m = nplotted_m + 1;
end

xlabel(axM,"\\mu (eV)","FontSize",FS);
ylabel(axM,sprintf("Re \\chi(%d,%d)", iq_pick, jq_pick),"FontSize",FS);
title(axM, sprintf("Re \\chi vs \\mu | T=%.6g K (tol %.1g) | curves=%d", ...
    T_target, T_tol, nplotted_m), "FontSize",FS, "FontWeight","normal");
grid(axM,"on"); box(axM,"on");
set(axM,"FontSize",FS,"LineWidth",1.0,"TickDir","out");

if use_colorbar
    caxis(axM,[Pmin Pmax]);
    cb = colorbar(axM);
    cb.Label.String = "polar\\_mu (meV)";
    cb.Label.Interpreter = "latex";
end

lg = legend(axM, "Location", "eastoutside");
lg.Title.String = "polar\\_mu";
lg.Title.Interpreter = "latex";
lg.FontSize = FS-2;

if save_png
    out.png_mu = fullfile(out_dir, sprintf("chi_vs_mu_T%.6g_i%d_j%d.png", T_target, iq_pick, jq_pick));
    exportgraphics(figM, out.png_mu, "Resolution", 220);
    fprintf("[saved] %s\n", out.png_mu);
end

% ============================================================
% nested helper: progress
% ============================================================
function progress_print()
    t_now = toc(t0);
    if (t_now - t_last) > print_every || k==Nraw
        pct = 100*k/Nraw;
        if k <= 1
            fprintf("  %5d/%5d (%.1f%%)  elapsed %.1fs\n", k, Nraw, pct, t_now);
        else
            eta = (t_now/k) * (Nraw-k);
            fprintf("  %5d/%5d (%.1f%%)  elapsed %.1fs  ETA %.1fs\n", k, Nraw, pct, t_now, eta);
        end
        drawnow;
        t_last = t_now;
    end
end

end % ===== end main =====


% ============================================================
% helpers (ALL included)
% ============================================================

function s = range_str(r)
if isempty(r)
    s = "[]";
else
    s = sprintf("[%.12g, %.12g]", min(r), max(r));
end
end

function tf = pass_range_AND(v, range)
if isempty(range)
    tf = true; return;
end
if ~isfinite(v)
    tf = false; return;
end
lo = min(range); hi = max(range);
tf = (v >= lo) && (v <= hi);
end

function polar = parse_polar_from_path(p)
% Accept:
%   .../polar0.01/...
%   .../polar_meV0.01/...
polar = NaN;
tok = regexp(char(p), "[\\/](?:polar_meV|polar)\s*([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)", "tokens", "once");
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

function add_point_by_polar(M, key, polar0, mu, dop, chi, path)
if ~isKey(M, key)
    S = struct();
    S.polar0 = polar0;
    S.mu  = mu;
    S.dop = dop;
    S.chi = chi;
    S.path = string(path);
    M(key) = S;
else
    S = M(key);
    S.mu(end+1,1)   = mu;
    S.dop(end+1,1)  = dop;
    S.chi(end+1,1)  = chi;
    S.path(end+1,1) = string(path);
    M(key) = S;
end
end

function [x2, y2] = collapse_same_x(x, y)
x = double(x(:)); y = double(y(:));
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
if isempty(x)
    x2 = []; y2 = []; return;
end
[xu, ~, ic] = unique(x);
ym = accumarray(ic, y, [], @mean);
[x2, ord] = sort(xu);
y2 = ym(ord);
end