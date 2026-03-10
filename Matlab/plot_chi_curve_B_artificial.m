function out = plot_chi_curve_B_artificial()
% plot_chi_curve_B_artificial_fromBaselineOnly
% Artificial chi(B) built ONLY from user-selected baseline polar0 folder.
% No reading of other polar folders.
%
% Steps:
%   1) UI choose baseline polar0 folder
%   2) read baseline chi0(mu) from that folder
%   3) user provides polar_list (meV)
%   4) chi_art(mu; p) = 0.5*(chi0(mu+p/2)+chi0(mu-p/2))
%   5) plot curves + plot max chi vs polar
%
% Output:
%   out.png_art_mu
%   out.png_max_vs_polar

out = struct();
out.png_art_mu = "";
out.png_max_vs_polar = "";

% ============================================================
% USER SETTINGS (EDIT HERE)
% ============================================================
D_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_200\D0.084";
% D_root = "/Users/haoranyan/Library/CloudStorage/OneDrive-Emory/Rhombohedral_SC/rhombohedral_project/data/chi_sk_mu_200/D0.067";

% ---- pick ONE temperature ----
T_target = 4;        % K
T_tol    = 1e-6;     % K tolerance

% ---- pick ONE (iq,jq) ----
iq_pick = -133;
jq_pick = -133;

% ---- baseline data filters (AND); [] => disable ----
mu_range    = [0, 1000];   % unit depends on mu_unit_in_header
dop_range   = [];          % [] => all doping

% ---- USER SPECIFIED polar list (meV) to construct artificial curves ----
% (include 0 if you want, it will reproduce baseline)
polar_list = [0, 2.2];

% ---- unit of mu in header ----
% "eV": polar_meV will be converted to eV when shifting (meV * 1e-3)
% "meV": shift directly
mu_unit_in_header = "eV";

% ---- output ----
save_png = true;
out_dir  = fullfile(string(D_root), "plot_Bcurves_artificial_baselineOnly");
if save_png && ~exist(out_dir,"dir"), mkdir(out_dir); end

% ---- style ----
FS = 14;
use_colorbar = false;
% ============================================================

assert(isfolder(D_root), "D_root not found: %s", D_root);

fprintf("Root D folder:\n  %s\n", D_root);
fprintf("Pick (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);
fprintf("Pick T_target=%.12g K, tol=%.3g K => range=[%.12g, %.12g]\n", ...
    T_target, T_tol, T_target-T_tol, T_target+T_tol);
fprintf("[baseline filter AND]\n");
fprintf("  mu_range    = %s (%s)\n", range_str(mu_range), mu_unit_in_header);
fprintf("  dop_range   = %s\n", range_str(dop_range));
fprintf("[artificial] polar_list (meV) = %s\n", mat2str(polar_list));

% ============================================================
% 1) UI: choose baseline polar0 folder
% ============================================================
baseline_dir = uigetdir(char(D_root), "Select BASELINE polar0 folder (use ONLY this folder to build chi_art)");
if isequal(baseline_dir, 0)
    error("User canceled baseline folder selection.");
end
baseline_dir = string(baseline_dir);
assert(isfolder(baseline_dir), "Selected baseline_dir not found: %s", baseline_dir);
fprintf("[baseline] selected folder:\n  %s\n", baseline_dir);

% Try to parse polar value from path (just for labeling)
polar0_from_path = parse_polar_from_path(baseline_dir);
if isfinite(polar0_from_path)
    fprintf("[baseline] polar parsed from path: %.12g meV\n", polar0_from_path);
else
    fprintf("[baseline] polar parsed from path: NaN (will just label as baseline)\n");
end

% ============================================================
% 2) Collect baseline curve chi0(mu) only from baseline folder
% ============================================================
[mu0, chi0] = collect_curve_from_folder(baseline_dir, T_target, T_tol, iq_pick, jq_pick, mu_range, dop_range);
if numel(mu0) < 3
    error("Baseline curve has too few points (need >=3). Check baseline folder / T_target / (iq,jq) / filters.");
end
fprintf("[baseline] points=%d, mu-range=[%.6g, %.6g] (%s)\n", numel(mu0), min(mu0), max(mu0), mu_unit_in_header);

mu0  = double(mu0(:));
chi0 = double(chi0(:));
F0   = griddedInterpolant(mu0, chi0, "linear", "none");

% unit conversion for shift
switch lower(string(mu_unit_in_header))
    case "ev"
        polar_to_mu = 1e-3;  % meV -> eV
    case "mev"
        polar_to_mu = 1.0;   % meV -> meV
    otherwise
        error('mu_unit_in_header must be "eV" or "meV".');
end

% sanitize polar_list
polar_list = unique(double(polar_list(:)));
polar_list = polar_list(isfinite(polar_list));
polar_list = sort(polar_list);

if isempty(polar_list)
    error("polar_list is empty.");
end

% ============================================================
% 3) Plot artificial chi(mu) curves + baseline
% ============================================================
cmap = turbo(256);
Pmin = min(polar_list);
Pmax = max(polar_list);
if Pmin == Pmax, Pmax = Pmin + 1e-12; end

set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
set(groot, "defaultTextInterpreter", "latex");

figA = figure("Color","w","Units","pixels","Position",[120 120 1200 750], ...
    "Name","artificial chi(mu) from baseline-only");
axA = axes(figA,"Position",[0.08 0.12 0.82 0.82]); hold(axA,"on");
colormap(axA, cmap);

% baseline dashed
plot(axA, mu0, chi0, "--", "LineWidth", 2.2, "DisplayName", "baseline (chosen folder)");

% store max chi vs polar
chi_max = nan(numel(polar_list),1);

n_art = 0;
for i = 1:numel(polar_list)
    pmev = polar_list(i);
    dmu = 0.5 * pmev * polar_to_mu;

    y1 = F0(mu0 + dmu);
    y2 = F0(mu0 - dmu);
    y  = 0.5 * (y1 + y2);

    ok = isfinite(y);
    if nnz(ok) < 3
        chi_max(i) = NaN;
        continue;
    end

    x = mu0(ok);
    y = y(ok);

    ci = color_index(pmev, Pmin, Pmax);
    plot(axA, x, y, "-", ...
        "LineWidth", 1.8, ...
        "Color", cmap(ci,:), ...
        "DisplayName", sprintf("artificial polar = %.6g\\,meV", pmev));
    n_art = n_art + 1;

    chi_max(i) = max(y);
end

xlabel(axA, sprintf("\\mu (%s)", mu_unit_in_header), "FontSize", FS);
ylabel(axA, sprintf("Re \\chi(%d,%d)", iq_pick, jq_pick), "FontSize", FS);
title(axA, sprintf("Artificial \\chi(\\mu) from baseline-only shift-average | T=%.6g K | curves=%d", ...
    T_target, n_art), "FontSize", FS, "FontWeight","normal");
grid(axA,"on"); box(axA,"on");
set(axA,"FontSize",FS,"LineWidth",1.0,"TickDir","out");

if use_colorbar
    caxis(axA,[Pmin Pmax]);
    cb = colorbar(axA);
    cb.Label.String = "polar\\_mu (meV)";
    cb.Label.Interpreter = "latex";
end

lg = legend(axA, "Location", "eastoutside");
lg.Title.String = "polar\\_mu";
lg.Title.Interpreter = "latex";
lg.FontSize = FS-2;

if save_png
    out.png_art_mu = fullfile(out_dir, sprintf("chi_artificial_vs_mu_T%.6g_i%d_j%d.png", T_target, iq_pick, jq_pick));
    exportgraphics(figA, out.png_art_mu, "Resolution", 220);
    fprintf("[saved] %s\n", out.png_art_mu);
end

% ============================================================
% 4) Plot max chi vs polar
% ============================================================
ok2 = isfinite(chi_max);
p2  = polar_list(ok2);
m2  = chi_max(ok2);

figB = figure("Color","w","Units","pixels","Position",[160 160 980 620], ...
    "Name","max artificial chi vs polar (baseline-only)");
axB = axes(figB,"Position",[0.10 0.14 0.82 0.80]); hold(axB,"on");

plot(axB, p2, m2, "o-", "LineWidth", 1.8, "MarkerSize", 6);

xlabel(axB, "polar\_mu (meV)", "FontSize", FS);
ylabel(axB, sprintf("\\max\\, Re\\chi(%d,%d)", iq_pick, jq_pick), "FontSize", FS);
title(axB, "Max of artificial \\chi vs polar (constructed from baseline only)", ...
    "FontSize", FS, "FontWeight","normal");
grid(axB,"on"); box(axB,"on");
set(axB,"FontSize",FS,"LineWidth",1.0,"TickDir","out");

if save_png
    out.png_max_vs_polar = fullfile(out_dir, sprintf("chi_artificial_max_vs_polar_T%.6g_i%d_j%d.png", T_target, iq_pick, jq_pick));
    exportgraphics(figB, out.png_max_vs_polar, "Resolution", 220);
    fprintf("[saved] %s\n", out.png_max_vs_polar);
end

end % ===== end main =====


% ============================================================
% helpers (ALL included)
% ============================================================

function [mu0, chi0] = collect_curve_from_folder(folder, T_target, T_tol, iq_pick, jq_pick, mu_range, dop_range)
files0 = dir(fullfile(folder, "**", "chi*.txt"));
files0 = files0(~[files0.isdir]);
if isempty(files0)
    error("No chi*.txt found under selected baseline folder: %s", folder);
end

T_lo = T_target - T_tol;
T_hi = T_target + T_tol;

mu_list = [];
chi_list = [];

for k = 1:numel(files0)
    p = string(fullfile(files0(k).folder, files0(k).name));
    m = parse_header_T_doping_mu(p);

    if ~isfinite(m.T), continue; end
    if ~(m.T >= T_lo && m.T <= T_hi), continue; end
    if ~pass_range_AND(m.mu, mu_range) || ~pass_range_AND(m.doping, dop_range), continue; end

    [iq_col, jq_col, re_col] = read_chi_table_cols(p);
    if isempty(iq_col), continue; end

    tf = (iq_col(:) == iq_pick) & (jq_col(:) == jq_pick);
    if ~any(tf), continue; end

    mu_list(end+1,1)  = m.mu; %#ok<AGROW>
    chi_list(end+1,1) = double(re_col(find(tf,1,'first'))); %#ok<AGROW>
end

[mu0, chi0] = collapse_same_x(mu_list, chi_list);
end

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

function ci = color_index(pval, Pmin, Pmax)
ci = round((pval - Pmin)/(Pmax - Pmin) * 255) + 1;
ci = max(min(ci,256),1);
end