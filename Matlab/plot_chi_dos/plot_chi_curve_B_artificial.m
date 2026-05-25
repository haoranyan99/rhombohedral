function out = plot_chi_curve_B_artificial()
% Artificial chi(B) built from selected baseline polar0 folder.
% Shift-average is done in mu space:
%   chi_art(mu; p) = 0.5 * [chi0(mu+p/2) + chi0(mu-p/2)]
% Plot x-axis uses doping parsed from the same baseline files.
%
% polar=0: gray dashed line
% polar>0: deep blue -> light blue as artificial polar increases

out = struct();
out.png_art_doping = "";
out.png_max_vs_polar = "";

% ============================================================
% USER SETTINGS
% ============================================================
D_root = "E:\rg_master\data";

T_target = 6.5;
T_tol    = 1e-6;

iq_pick = 3;
jq_pick = -3;

mu_range  = [];
dop_range = [-2.5, 0];

polar_list = 0:1:6;          % meV
mu_unit_in_header = "eV";    % "eV" or "meV"

save_png = true;
out_dir  = fullfile(string(D_root), "plot_Bcurves_artificial_baselineOnly");
if save_png && ~exist(out_dir,"dir")
    mkdir(out_dir);
end

FS = 14;
use_colorbar = true;

% ============================================================
% Select baseline folder
% ============================================================
assert(isfolder(D_root), "D_root not found: %s", D_root);

baseline_dir = uigetdir(char(D_root), ...
    "Select BASELINE polar0 folder");

if isequal(baseline_dir, 0)
    error("User canceled baseline folder selection.");
end

baseline_dir = string(baseline_dir);

% ============================================================
% Collect baseline curve: mu, doping, chi
% ============================================================
[mu0, dop0, chi0] = collect_curve_from_folder( ...
    baseline_dir, T_target, T_tol, iq_pick, jq_pick, mu_range, dop_range);

if numel(mu0) < 3
    error("Baseline curve has too few points.");
end

mu0  = double(mu0(:));
dop0 = double(dop0(:));
chi0 = double(chi0(:));

% sort by mu
[mu0, ord] = sort(mu0);
dop0 = dop0(ord);
chi0 = chi0(ord);

F0 = griddedInterpolant(mu0, chi0, "linear", "none");

switch lower(string(mu_unit_in_header))
    case "ev"
        polar_to_mu = 1e-3;
    case "mev"
        polar_to_mu = 1.0;
    otherwise
        error('mu_unit_in_header must be "eV" or "meV".');
end

polar_list = unique(double(polar_list(:)));
polar_list = polar_list(isfinite(polar_list));
polar_list = sort(polar_list);

if isempty(polar_list)
    error("polar_list is empty.");
end

Pmin = min(polar_list);
Pmax = max(polar_list);
if Pmin == Pmax
    Pmax = Pmin + 1e-12;
end

cmap = blue_gradient(256);

set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
set(groot, "defaultTextInterpreter", "latex");

% ============================================================
% Plot artificial chi curves vs doping
% ============================================================
figA = figure("Color","w", ...
    "Units","pixels", ...
    "Position",[120 120 1100 720], ...
    "Name","artificial chi(doping)");

axA = axes(figA, "Position",[0.10 0.13 0.78 0.80]);
hold(axA, "on");

chi_max = nan(numel(polar_list),1);

for i = 1:numel(polar_list)

    pmev = polar_list(i);
    dmu = 0.5 * pmev * polar_to_mu;

    y1 = F0(mu0 + dmu);
    y2 = F0(mu0 - dmu);
    y  = 0.5 * (y1 + y2);

    ok = isfinite(y) & isfinite(dop0);

    if nnz(ok) < 3
        continue;
    end

    x = dop0(ok);
    y = y(ok);

    chi_max(i) = max(y);

    if abs(pmev) < 1e-12
        plot(axA, x, y, "--", ...
            "Color", [0.55 0.55 0.55], ...
            "LineWidth", 2.4, ...
            "DisplayName", "$B_{\rm art}=0$");
    else
        ci = color_index(pmev, Pmin, Pmax);
        plot(axA, x, y, "-", ...
            "Color", cmap(ci,:), ...
            "LineWidth", 1.9, ...
            "HandleVisibility","off");
    end
end

xlabel(axA, "Doping $(10^{12}\,\mathrm{cm}^{-2})$", "FontSize", FS);
ylabel(axA, sprintf("$\\mathrm{Re}\\,\\chi(%d,%d)$", iq_pick, jq_pick), ...
    "FontSize", FS);

title(axA, sprintf("$T=%.3f\\,\\mathrm{K}$", T_target), ...
    "FontSize", FS, "FontWeight","normal");

box(axA, "on");
grid(axA, "off");

set(axA, ...
    "FontSize", FS, ...
    "LineWidth", 1.3, ...
    "TickDir", "in");

if use_colorbar
    colormap(axA, cmap);
    caxis(axA, [Pmin Pmax]);

    cb = colorbar(axA);
    cb.Label.String = "$B_{\rm art}$ / polar (meV)";
    cb.Label.Interpreter = "latex";
    cb.TickLabelInterpreter = "latex";
    cb.TickDirection = "in";
    cb.LineWidth = 1.2;
end

legend(axA, "Location", "best", "Box", "off");

if save_png
    out.png_art_doping = fullfile(out_dir, ...
        sprintf("chi_artificial_blue_vs_doping_T%.6g_i%d_j%d.png", ...
        T_target, iq_pick, jq_pick));
    exportgraphics(figA, out.png_art_doping, "Resolution", 220);
end

% ============================================================
% Plot max chi vs polar
% ============================================================
ok2 = isfinite(chi_max);
p2  = polar_list(ok2);
m2  = chi_max(ok2);

figB = figure("Color","w", ...
    "Units","pixels", ...
    "Position",[160 160 900 620], ...
    "Name","max chi vs artificial polar");

axB = axes(figB, "Position",[0.12 0.15 0.80 0.78]);
hold(axB, "on");

plot(axB, p2, m2, "o-", ...
    "Color", [0.05 0.20 0.60], ...
    "LineWidth", 1.9, ...
    "MarkerSize", 6, ...
    "MarkerFaceColor", [0.70 0.85 1.00]);

xlabel(axB, "$B_{\rm art}$ / polar (meV)", "FontSize", FS);
ylabel(axB, sprintf("$\\max\\,\\mathrm{Re}\\,\\chi(%d,%d)$", ...
    iq_pick, jq_pick), "FontSize", FS);

title(axB, sprintf("$T=%.3f\\,\\mathrm{K}$", T_target), ...
    "FontSize", FS, "FontWeight","normal");

box(axB, "on");
grid(axB, "off");

set(axB, ...
    "FontSize", FS, ...
    "LineWidth", 1.3, ...
    "TickDir", "in");

if save_png
    out.png_max_vs_polar = fullfile(out_dir, ...
        sprintf("chi_artificial_max_vs_polar_T%.6g_i%d_j%d.png", ...
        T_target, iq_pick, jq_pick));
    exportgraphics(figB, out.png_max_vs_polar, "Resolution", 220);
end

fprintf("\nSaved curve figure:\n%s\n", out.png_art_doping);
fprintf("Saved max figure:\n%s\n", out.png_max_vs_polar);

end

% ============================================================
% Helpers
% ============================================================

function [mu0, dop0, chi0] = collect_curve_from_folder(folder, T_target, T_tol, ...
    iq_pick, jq_pick, mu_range, dop_range)

files0 = dir(fullfile(folder, "**", "chi*.txt"));
files0 = files0(~[files0.isdir]);

if isempty(files0)
    error("No chi*.txt found under selected folder: %s", folder);
end

T_lo = T_target - T_tol;
T_hi = T_target + T_tol;

mu_list  = [];
dop_list = [];
chi_list = [];

for k = 1:numel(files0)

    p = string(fullfile(files0(k).folder, files0(k).name));
    m = parse_header_T_doping_mu(p);

    if ~isfinite(m.T), continue; end
    if ~(m.T >= T_lo && m.T <= T_hi), continue; end
    if ~pass_range_AND(m.mu, mu_range), continue; end
    if ~pass_range_AND(m.doping, dop_range), continue; end

    [iq_col, jq_col, re_col] = read_chi_table_cols(p);
    if isempty(iq_col), continue; end

    tf = (iq_col(:) == iq_pick) & (jq_col(:) == jq_pick);
    if ~any(tf), continue; end

    mu_list(end+1,1)  = m.mu; %#ok<AGROW>
    dop_list(end+1,1) = m.doping; %#ok<AGROW>
    chi_list(end+1,1) = double(re_col(find(tf,1,'first'))); %#ok<AGROW>
end

[mu0, dop0, chi0] = collapse_same_mu(mu_list, dop_list, chi_list);

end

function [mu2, dop2, chi2] = collapse_same_mu(mu, dop, chi)

mu  = double(mu(:));
dop = double(dop(:));
chi = double(chi(:));

ok = isfinite(mu) & isfinite(dop) & isfinite(chi);
mu  = mu(ok);
dop = dop(ok);
chi = chi(ok);

if isempty(mu)
    mu2 = [];
    dop2 = [];
    chi2 = [];
    return;
end

[mu_unique, ~, ic] = unique(mu);

dop_mean = accumarray(ic, dop, [], @mean);
chi_mean = accumarray(ic, chi, [], @mean);

[mu2, ord] = sort(mu_unique);
dop2 = dop_mean(ord);
chi2 = chi_mean(ord);

end

function tf = pass_range_AND(v, range)

if isempty(range)
    tf = true;
    return;
end

if ~isfinite(v)
    tf = false;
    return;
end

lo = min(range);
hi = max(range);
tf = (v >= lo) && (v <= hi);

end

function m = parse_header_T_doping_mu(in_path)

m = struct("T",NaN, "doping",NaN, "mu",NaN);

fid = fopen(in_path,'r');
if fid < 0
    return;
end

cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

while true
    tline = fgetl(fid);
    if ~ischar(tline), break; end

    s0 = strtrim(tline);
    if ~startsWith(s0, "#"), break; end

    s = strtrim(erase(s0, "#"));

    if ~isfinite(m.T)
        tok = regexp(s, "(?:^|\s)T(?:_K)?\s*=\s*" + num, "tokens", "once");
        if ~isempty(tok)
            m.T = str2double(tok{1});
        end
    end

    if ~isfinite(m.doping)
        tok = regexp(s, "doping\s*=\s*" + num, "tokens", "once");
        if ~isempty(tok)
            m.doping = str2double(tok{1});
        end
    end

    if ~isfinite(m.mu)
        tok = regexp(s, "(?:^|\s)(?:mu|mu_eV|EF|EF_eV)\s*=\s*" + num, ...
            "tokens", "once");
        if ~isempty(tok)
            m.mu = str2double(tok{1});
        end
    end
end

end

function [iq, jq, re] = read_chi_table_cols(in_path)

iq = [];
jq = [];
re = [];

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

function ci = color_index(pval, Pmin, Pmax)

ci = round((pval - Pmin) / (Pmax - Pmin) * 255) + 1;
ci = max(min(ci, 256), 1);

end

function cmap = blue_gradient(N)

if nargin < 1
    N = 256;
end

% 🔥 更极端的深蓝和浅蓝（增强对比）
dark_blue  = [0.00, 0.05, 0.35];   % 非常深的蓝（接近navy）
light_blue = [0.85, 0.93, 1.00];   % 非常浅的蓝（接近白）

t = linspace(0,1,N)';

% 🔥 非线性拉伸（关键！）
gamma = 0.6;   % <1：增强浅色变化
t = t.^gamma;

cmap = [ ...
    dark_blue(1)  + (light_blue(1)-dark_blue(1)) * t, ...
    dark_blue(2)  + (light_blue(2)-dark_blue(2)) * t, ...
    dark_blue(3)  + (light_blue(3)-dark_blue(3)) * t ...
];

end