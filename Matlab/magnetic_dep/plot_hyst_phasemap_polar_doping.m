function plot_hyst_phasemap_polar_doping()
clc; close all;

% =========================
% USER SETTINGS
% =========================
T_target = 6.5;          % fixed temperature
quantity = "psi";        % "psi" or "X"
FS = 16;

doping_N = 300;          % common doping grid number
use_abs_diff = true;     % true: |fwd-bwd|, false: fwd-bwd

default_dir = 'E:\rg_master\rhombohedral\Matlab\hyst_data';

% =========================
% Select MAT file
% =========================
if ~isfolder(default_dir)
    default_dir = pwd;
end

[fname, fpath] = uigetfile(fullfile(default_dir, '*.mat'), ...
    'Select hyst_allT_allArtificialPolar.mat');

if isequal(fname,0)
    return;
end

matfile = fullfile(fpath, fname);
S = load(matfile);

if ~isfield(S, "HYST")
    error("Selected MAT file does not contain variable HYST.");
end

HYST = S.HYST;

valid = arrayfun(@(x) isfield(x,"T") && isfinite(x.T) && ...
    isfield(x,"artificial_polar") && isfinite(x.artificial_polar), HYST);

entries = HYST(valid);

if isempty(entries)
    error("No valid HYST entries found.");
end

% =========================
% Select fixed temperature
% =========================
tolT = 1e-8;
E = entries(abs([entries.T] - T_target) < tolT);

if isempty(E)
    Tall = unique([entries.T]);
    error("No data found for T = %.6g K. Available T: %s", ...
        T_target, mat2str(Tall));
end

% =========================
% Build common doping grid
% =========================
all_dop_min = -inf;
all_dop_max = inf;

for i = 1:numel(E)
    [xf, ~, xb, ~] = get_quantity_curves_(E(i), quantity);

    xf = xf(isfinite(xf));
    xb = xb(isfinite(xb));

    if isempty(xf) || isempty(xb)
        continue;
    end

    all_dop_min = max(all_dop_min, max(min(xf), min(xb)));
    all_dop_max = min(all_dop_max, min(max(xf), max(xb)));
end

if ~isfinite(all_dop_min) || ~isfinite(all_dop_max) || all_dop_min >= all_dop_max
    error("Cannot build common doping range.");
end

dop_grid = linspace(all_dop_min, all_dop_max, doping_N);

% =========================
% Compute phase map
% =========================
polar = nan(numel(E),1);
Z = nan(numel(E), doping_N);

for i = 1:numel(E)
    R = E(i);
    polar(i) = R.artificial_polar;

    [xf, yf, xb, yb] = get_quantity_curves_(R, quantity);

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

okrow = isfinite(polar) & any(isfinite(Z),2);
polar = polar(okrow);
Z = Z(okrow,:);

[polar, ord] = sort(polar);
Z = Z(ord,:);

if isempty(polar)
    error("No valid polar-doping phase map data found.");
end

% =========================
% Plot
% =========================
figure('Color','w','Position',[120 120 850 620]);

imagesc(dop_grid, polar, Z);
set(gca, 'YDir', 'normal');

colormap(flipud(hot));
cb = colorbar;
cb.TickLabelInterpreter = 'latex';

if use_abs_diff
    cb.Label.String = sprintf('$||%s|_{\\rm fwd} - |%s|_{\\rm bwd}|$', quantity, quantity);
else
    cb.Label.String = sprintf('$|%s|_{\\rm fwd} - |%s|_{\\rm bwd}$', quantity, quantity);
end
cb.Label.Interpreter = 'latex';

xlabel('Doping $(10^{12}\,\mathrm{cm}^{-2})$', 'Interpreter','latex');
ylabel('Artificial polar / magnetic field proxy (meV)', 'Interpreter','latex');

title(sprintf('$T = %.3f\\,\\mathrm{K}$, hysteresis phase map of $|%s|$', ...
    T_target, quantity), 'Interpreter','latex');

set(gca, ...
    'FontSize', FS, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth', 1.3, ...
    'TickDir','out', ...
    'Box','on');

fprintf("\nLoaded: %s\n", matfile);
fprintf("T = %.6g K\n", T_target);
fprintf("Quantity = %s\n", quantity);
fprintf("Doping range = [%.6g, %.6g]\n", dop_grid(1), dop_grid(end));
fprintf("Polar range = [%.6g, %.6g]\n", min(polar), max(polar));

end

function [xf, yf, xb, yb] = get_quantity_curves_(R, quantity)

xf = R.dop_fwd(:);
xb = R.dop_bwd(:);

switch lower(string(quantity))
    case "psi"
        if isfield(R,"psi_f_plot") && isfield(R,"psi_b_plot")
            yf = R.psi_f_plot(:);
            yb = R.psi_b_plot(:);
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
            yf = R.X_f_plot(:);
            yb = R.X_b_plot(:);
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