function plot_hyst_phasemap_polar_doping()

clc; close all;

% =========================
% USER SETTINGS
% =========================
T_target = 6.5;
quantity = "psi";

FS = 16;
doping_N = 300;
use_abs_diff = true;

white_pos = 0.25;   % 🔥 手动控制白色位置（推荐 0.02~0.1）

default_dir = 'E:\rg_master\rhombohedral\Matlab\hyst_data';

% =========================
% Load MAT
% =========================
if ~isfolder(default_dir)
    default_dir = pwd;
end

[fname, fpath] = uigetfile(fullfile(default_dir, '*.mat'), ...
    'Select hyst_allT_allArtificialPolar.mat');

if isequal(fname,0); return; end

S = load(fullfile(fpath, fname));
HYST = S.HYST;

valid = arrayfun(@(x) isfield(x,"T") && isfield(x,"artificial_polar"), HYST);
entries = HYST(valid);

tolT = 1e-8;
E = entries(abs([entries.T] - T_target) < tolT);

if isempty(E)
    error("No data at T = %.6g", T_target);
end

% =========================
% Build doping grid
% =========================
dop_min = -inf;
dop_max = inf;

for i = 1:numel(E)
    [xf,~,xb,~] = get_quantity(E(i), quantity);

    xf = xf(isfinite(xf));
    xb = xb(isfinite(xb));

    if isempty(xf) || isempty(xb); continue; end

    dop_min = max(dop_min, max(min(xf), min(xb)));
    dop_max = min(dop_max, min(max(xf), max(xb)));
end

dop_grid = linspace(dop_min, dop_max, doping_N);

% =========================
% Compute Z
% =========================
polar = nan(numel(E),1);
Z = nan(numel(E), doping_N);

for i = 1:numel(E)

    R = E(i);
    polar(i) = R.artificial_polar;

    [xf,yf,xb,yb] = get_quantity(R, quantity);

    okf = isfinite(xf) & isfinite(yf);
    okb = isfinite(xb) & isfinite(yb);

    xf = xf(okf); yf = yf(okf);
    xb = xb(okb); yb = yb(okb);

    if numel(xf) < 3 || numel(xb) < 3; continue; end

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
% sort
% =========================
[polar, ord] = sort(polar);
Z = Z(ord,:);

% =========================
% Plot
% =========================
figure('Color','w','Position',[120 120 850 620]);

imagesc(dop_grid, polar, Z);
set(gca, 'YDir','normal');

% 🔥 RWB colormap
colormap(red_white_blue_skewed(256, white_pos));

% 🔥 caxis（关键）
if use_abs_diff
    zmin = 0;
else
    zmin = min(Z(:), [], 'omitnan');
end
zmax = max(Z(:), [], 'omitnan');

if ~(zmax > zmin)
    zmax = zmin + eps;
end

caxis([zmin, zmax]);

cb = colorbar;
cb.Label.String = sprintf('$||%s|_{fwd}-|%s|_{bwd}|$', quantity, quantity);
cb.Label.Interpreter = 'latex';

xlabel('Doping', 'Interpreter','latex');
ylabel('Polar (meV)', 'Interpreter','latex');

title(sprintf('$T = %.2f\\,K$', T_target), 'Interpreter','latex');

set(gca, ...
    'FontSize', FS, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth', 1.3, ...
    'TickDir','out', ...
    'Box','on');

end

% =========================
% helper: get quantity
% =========================
function [xf,yf,xb,yb] = get_quantity(R, quantity)

xf = R.dop_fwd(:);
xb = R.dop_bwd(:);

switch lower(string(quantity))
    case "psi"
        yf = R.psi_f_plot(:);
        yb = R.psi_b_plot(:);
    case "x"
        yf = R.X_f_plot(:);
        yb = R.X_b_plot(:);
    otherwise
        error('quantity must be psi or X');
end

end

% =========================
% 🔥 RED-WHITE-BLUE
% =========================
function cmap = red_white_blue_skewed(N, white_pos)

if nargin < 1, N = 256; end
if nargin < 2, white_pos = 0.05; end

blue  = [45, 75, 145] / 255;
white = [1, 1, 1];
red   = [235, 70, 50] / 255;

n_blue = max(2, round(N * white_pos));
n_red  = N - n_blue;

c1 = [linspace(blue(1),  white(1), n_blue)', ...
      linspace(blue(2),  white(2), n_blue)', ...
      linspace(blue(3),  white(3), n_blue)'];

c2 = [linspace(white(1), red(1), n_red)', ...
      linspace(white(2), red(2), n_red)', ...
      linspace(white(3), red(3), n_red)'];

cmap = [c1; c2];

end