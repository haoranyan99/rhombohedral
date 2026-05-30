function OUT = plot_hyst_strength_map_B_doping()
clc; close all;

% =========================
% USER SETTINGS
% =========================
T_target = 6.4;
artificial_polar_target = 0;
quantity = "psi";          % "psi" or "X"

doping_range = [-1, 0];    % 10^12 cm^-2
doping_N = 300;

use_abs_diff = true;       % true: abs(fwd-bwd)
FS = 16;

default_dir = 'E:\rg_master\rhombohedral\Matlab\hyst';
if ~isfolder(default_dir)
    default_dir = 'E:\rg_master\rhombohedral\Matlab\hyst_data';
end

% =========================
% Select MAT file
% =========================
if ~isfolder(default_dir)
    default_dir = pwd;
end

[fname, fpath] = uigetfile(fullfile(default_dir, 'hyst*.mat'), ...
    'Select hyst_allB*.mat');

if isequal(fname,0)
    OUT = [];
    return;
end

matfile = fullfile(fpath, fname);
S = load(matfile);

if ~isfield(S, "HYST")
    error("Selected MAT file does not contain variable HYST.");
end

HYST = S.HYST;

% =========================
% Select valid entries
% =========================
valid = arrayfun(@is_valid_hyst_entry_, HYST);
entries = HYST(valid);

if isempty(entries)
    error("No valid real-B HYST entries found.");
end

tolT = 1e-8;
tolP = 1e-8;

E = entries(abs([entries.T] - T_target) < tolT);
if isempty(E)
    Tall = unique([entries.T]);
    error("No data found for T = %.6g K. Available T: %s", ...
        T_target, mat2str(Tall));
end

if isfield(E, "artificial_polar")
    E = E(abs([E.artificial_polar] - artificial_polar_target) < tolP);
    if isempty(E)
        error("No data found for artificial polar = %.6g meV.", ...
            artificial_polar_target);
    end
end

% =========================
% Common doping grid
% =========================
dop_grid = linspace(doping_range(1), doping_range(2), doping_N);

B_list = nan(numel(E),1);
Z = nan(numel(E), doping_N);

% =========================
% Compute hysteresis strength map
% =========================
for i = 1:numel(E)

    R = E(i);
    B_list(i) = R.B_T;

    [xf, yf, xb, yb] = get_quantity_curves_(R, quantity);

    xf = xf(:); yf = yf(:);
    xb = xb(:); yb = yb(:);

    okf = isfinite(xf) & isfinite(yf);
    okb = isfinite(xb) & isfinite(yb);

    xf = xf(okf); yf = yf(okf);
    xb = xb(okb); yb = yb(okb);

    if numel(xf) < 3 || numel(xb) < 3
        continue;
    end

    [xf, ord] = sort(xf); yf = yf(ord);
    [xb, ord] = sort(xb); yb = yb(ord);

    [xf, ia] = unique(xf); yf = yf(ia);
    [xb, ia] = unique(xb); yb = yb(ia);

    yf_grid = interp1(xf, yf, dop_grid, 'linear', NaN);
    yb_grid = interp1(xb, yb, dop_grid, 'linear', NaN);

    if use_abs_diff
        Z(i,:) = abs(yf_grid - yb_grid);
    else
        Z(i,:) = yf_grid - yb_grid;
    end
end

% =========================
% Clean and sort
% =========================
okrow = isfinite(B_list) & any(isfinite(Z),2);

B_list = B_list(okrow);
Z = Z(okrow,:);

[B_list, ord] = sort(B_list);
Z = Z(ord,:);

if isempty(B_list)
    error("No valid B-doping hysteresis map data found.");
end

% =========================
% Plot
% =========================
fig = figure('Color','w','Position',[120 120 880 620]);
ax = axes(fig);

imagesc(ax, dop_grid, B_list, Z);
set(ax, 'YDir', 'normal');

colormap(ax, hot);
cb = colorbar(ax);
cb.TickLabelInterpreter = 'latex';
cb.LineWidth = 1.2;

if use_abs_diff
    cb.Label.String = sprintf('$\\Delta |%s| = ||%s|_{\\rm fwd}-|%s|_{\\rm bwd}|$', ...
        quantity, quantity, quantity);
else
    cb.Label.String = sprintf('$\\Delta |%s| = |%s|_{\\rm fwd}-|%s|_{\\rm bwd}$', ...
        quantity, quantity, quantity);
end
cb.Label.Interpreter = 'latex';

xlabel(ax, 'Doping $(10^{12}\,\mathrm{cm}^{-2})$', 'Interpreter','latex');
ylabel(ax, '$B$ (T)', 'Interpreter','latex');

title(ax, sprintf('$T=%.3f\\,\\mathrm{K}$, hysteresis strength map of $|%s|$', ...
    T_target, quantity), ...
    'Interpreter','latex', ...
    'FontWeight','normal');

set(ax, ...
    'FontSize', FS, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth', 1.3, ...
    'TickDir','out', ...
    'Box','on');

% =========================
% Output
% =========================
OUT = struct();
OUT.matfile = matfile;
OUT.T_target = T_target;
OUT.artificial_polar_target = artificial_polar_target;
OUT.quantity = quantity;
OUT.B_list = B_list;
OUT.dop_grid = dop_grid;
OUT.Z = Z;
OUT.fig = fig;
OUT.ax = ax;

fprintf("\nLoaded: %s\n", matfile);
fprintf("T = %.6g K\n", T_target);
fprintf("Artificial polar = %.6g meV\n", artificial_polar_target);
fprintf("Quantity = %s\n", quantity);
fprintf("Doping range = [%.6g, %.6g]\n", dop_grid(1), dop_grid(end));
fprintf("B range = [%.6g, %.6g] T\n", min(B_list), max(B_list));

end

% ============================================================
% Helpers
% ============================================================
function ok = is_valid_hyst_entry_(R)

ok = isfield(R,"T") && isfield(R,"B_T") && ...
    isfinite(R.T) && isfinite(R.B_T) && ...
    isfield(R,"dop_fwd") && isfield(R,"dop_bwd");

end

function [xf, yf, xb, yb] = get_quantity_curves_(R, quantity)

xf = R.dop_fwd(:);
xb = R.dop_bwd(:);

switch lower(string(quantity))

    case "psi"
        if isfield(R,"psi_f_plot") && isfield(R,"psi_b_plot")
            yf = abs(R.psi_f_plot(:));
            yb = abs(R.psi_b_plot(:));
        elseif isfield(R,"abspsi_fwd") && isfield(R,"abspsi_bwd")
            yf = R.abspsi_fwd(:);
            yb = R.abspsi_bwd(:);
        elseif isfield(R,"psi_fwd") && isfield(R,"psi_bwd")
            yf = abs(R.psi_fwd(:));
            yb = abs(R.psi_bwd(:));
        else
            error("Cannot find psi data in HYST entry.");
        end

    case "x"
        if isfield(R,"X_f_plot") && isfield(R,"X_b_plot")
            yf = abs(R.X_f_plot(:));
            yb = abs(R.X_b_plot(:));
        elseif isfield(R,"absX_fwd") && isfield(R,"absX_bwd")
            yf = R.absX_fwd(:);
            yb = R.absX_bwd(:);
        elseif isfield(R,"X_fwd") && isfield(R,"X_bwd")
            yf = abs(R.X_fwd(:));
            yb = abs(R.X_bwd(:));
        else
            error("Cannot find X data in HYST entry.");
        end

    otherwise
        error("quantity must be psi or X.");
end

end
