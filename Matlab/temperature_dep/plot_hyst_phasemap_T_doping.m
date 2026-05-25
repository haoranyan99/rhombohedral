function plot_hyst_phasemap_T_doping()

clc; close all;

% =========================
% USER SETTINGS
% =========================
polar_target = 0.0;     % meV

quantity = "psi";

FS = 16;
doping_N = 300;

use_abs_diff = true;

hot_yellow_value_frac = 0.5;
hot_yellow_color_pos = 0.5;
hot_white_value = 1.2;
hot_white_color_pos = 0.06;
hot_red_value = 3.9;
hot_red_color_pos = 0.65;
white_pos = 0.75;       % used only for signed fwd-bwd maps

default_dir = 'E:\rg_master\rhombohedral\Matlab\hyst_data';

% =========================
% Load MAT
% =========================
if ~isfolder(default_dir)
    default_dir = pwd;
end

[fname, fpath] = uigetfile(fullfile(default_dir, '*.mat'), ...
    'Select hyst_allT_allArtificialPolar.mat');

if isequal(fname,0)
    return;
end

S = load(fullfile(fpath, fname));
HYST = S.HYST;

valid = arrayfun(@(x) ...
    isfield(x,"T") && isfield(x,"artificial_polar"), ...
    HYST);

entries = HYST(valid);

% =========================
% Select fixed polar
% =========================
tolP = 1e-8;

E = entries(abs([entries.artificial_polar] - polar_target) < tolP);

if isempty(E)
    error("No data at polar = %.6g meV", polar_target);
end

% =========================
% Build common doping grid
% =========================
dop_min = -inf;
dop_max = inf;

for i = 1:numel(E)

    [xf,~,xb,~] = get_quantity(E(i), quantity);

    xf = xf(isfinite(xf));
    xb = xb(isfinite(xb));

    if isempty(xf) || isempty(xb)
        continue;
    end

    dop_min = max(dop_min, max(min(xf), min(xb)));
    dop_max = min(dop_max, min(max(xf), max(xb)));

end

dop_grid = linspace(dop_min, dop_max, doping_N);

% =========================
% Compute Z(T,doping)
% =========================
Tlist = nan(numel(E),1);
Z = nan(numel(E), doping_N);

for i = 1:numel(E)

    R = E(i);

    Tlist(i) = R.T;

    [xf,yf,xb,yb] = get_quantity(R, quantity);

    okf = isfinite(xf) & isfinite(yf);
    okb = isfinite(xb) & isfinite(yb);

    xf = xf(okf);
    yf = yf(okf);

    xb = xb(okb);
    yb = yb(okb);

    if numel(xf) < 3 || numel(xb) < 3
        continue;
    end

    [xf, ia] = unique(xf);
    yf = yf(ia);

    [xb, ia] = unique(xb);
    yb = yb(ia);

    yf_grid = interp1(xf, yf, dop_grid, 'linear', NaN);
    yb_grid = interp1(xb, yb, dop_grid, 'linear', NaN);

    if use_abs_diff
        Z(i,:) = abs(yf_grid - yb_grid);
    else
        Z(i,:) = yf_grid - yb_grid;
    end

end

% =========================
% sort by temperature
% =========================
[Tlist, ord] = sort(Tlist);

Z = Z(ord,:);

% =========================
% Plot
% =========================
figure('Color','w','Position',[120 120 850 620]);

imagesc(dop_grid, Tlist, Z);

set(gca, 'YDir','normal');

if use_abs_diff
    zmin = 0;
else
    zmin = min(Z(:), [], 'omitnan');
end

zmax = max(Z(:), [], 'omitnan');

if ~(zmax > zmin)
    zmax = zmin + eps;
end

if use_abs_diff
    zred = min(max(hot_red_value, zmin), zmax);
    zyellow_auto = zmin + hot_yellow_value_frac * (zmax - zmin);
    zyellow = min(zyellow_auto, hot_white_value + 0.5 * (zred - hot_white_value));
    zyellow = min(max(zyellow, zmin), zmax);
    white_value_frac = (hot_white_value - zmin) / (zmax - zmin);
    yellow_value_frac = (zyellow - zmin) / (zmax - zmin);
    red_value_frac = (zred - zmin) / (zmax - zmin);
    positive_colormap = positive_hot_colormap( ...
        256, white_value_frac, hot_white_color_pos, ...
        yellow_value_frac, hot_yellow_color_pos, ...
        red_value_frac, hot_red_color_pos);
    clim([zmin zmax]);
else
    zyellow = NaN;
    clim([zmin zmax]);
end

cb = colorbar;

if use_abs_diff
    colormap(gca, positive_colormap);
    [tick_pos, tick_label] = positive_hot_colorbar_ticks( ...
        zmin, zyellow, zmax);
    cb.Ticks = tick_pos;
    cb.TickLabels = tick_label;
else
    colormap(gca, red_white_blue_skewed(256, white_pos));
end

cb.Label.String = ...
    sprintf('$||%s|_{fwd}-|%s|_{bwd}|$', ...
    quantity, quantity);

cb.Label.Interpreter = 'latex';

xlabel('Doping', 'Interpreter','latex');

ylabel('$T$ (K)', 'Interpreter','latex');

title(sprintf('Polar = %.3f meV', polar_target), ...
    'Interpreter','latex');

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
% RED-WHITE-BLUE
% =========================
function cmap = red_white_blue_skewed(N, white_pos)

if nargin < 1
    N = 256;
end
if nargin < 2
    white_pos = 0.05;
end

blue  = [45, 75, 145] / 255;
white = [1, 1, 1];
red   = [235, 70, 50] / 255;

n_blue = max(2, round(N * white_pos));
n_red  = N - n_blue;

c1 = [ ...
    linspace(blue(1),  white(1), n_blue)', ...
    linspace(blue(2),  white(2), n_blue)', ...
    linspace(blue(3),  white(3), n_blue)' ...
    ];

c2 = [ ...
    linspace(white(1), red(1), n_red)', ...
    linspace(white(2), red(2), n_red)', ...
    linspace(white(3), red(3), n_red)' ...
    ];

cmap = [c1; c2];

end

function cmap = positive_hot_colormap(N, white_value_frac, white_color_pos, yellow_value_frac, yellow_color_pos, red_value_frac, red_color_pos)

if nargin < 1
    N = 256;
end
if nargin < 2
    white_value_frac = 0.25;
end
if nargin < 3
    white_color_pos = 0.06;
end
if nargin < 4
    yellow_value_frac = 0.75;
end
if nargin < 5
    yellow_color_pos = 0.25;
end
if nargin < 6
    red_value_frac = 0.9;
end
if nargin < 7
    red_color_pos = 0.65;
end

white_value_frac = min(max(white_value_frac, 1e-6), 1 - 2e-6);
white_color_pos = min(max(white_color_pos, 1e-6), 1 - 2e-6);
yellow_value_frac = min(max(yellow_value_frac, 1e-6), 1 - 1e-6);
yellow_color_pos = min(max(yellow_color_pos, white_color_pos + 1e-6), 1 - 1e-6);
if yellow_value_frac <= white_value_frac
    yellow_value_frac = min(1 - 1e-6, white_value_frac + 1e-6);
end
red_value_frac = min(max(red_value_frac, yellow_value_frac + 1e-6), 1 - 1e-6);
red_color_pos = min(max(red_color_pos, yellow_color_pos + 1e-6), 1 - 1e-6);

base = flipud(hot(max(1024, 4 * N)));
r = linspace(0, 1, N).';
u = nan(size(r));

lo = r <= white_value_frac;
mid1 = r > white_value_frac & r <= yellow_value_frac;
mid2 = r > yellow_value_frac & r <= red_value_frac;
hi = r > red_value_frac;
u(lo) = white_color_pos * r(lo) / white_value_frac;
u(mid1) = white_color_pos + (yellow_color_pos - white_color_pos) * ...
    (r(mid1) - white_value_frac) / (yellow_value_frac - white_value_frac);
u(mid2) = yellow_color_pos + (red_color_pos - yellow_color_pos) * ...
    (r(mid2) - yellow_value_frac) / (red_value_frac - yellow_value_frac);
u(hi) = red_color_pos + (1 - red_color_pos) * ...
    (r(hi) - red_value_frac) / (1 - red_value_frac);

xbase = linspace(0, 1, size(base, 1)).';
cmap = interp1(xbase, base, u, "linear");

end

function [ticks, labels] = positive_hot_colorbar_ticks(zmin, zyellow, zmax)

ticks = unique([ceil(zmin), round(zyellow), floor(zmax)]);
ticks = ticks(ticks >= zmin & ticks <= zmax);

if numel(ticks) < 3
    ticks = unique([zmin, zyellow, zmax]);
    labels = cellstr(compose('%.4g', ticks));
else
    labels = cellstr(compose('%.0f', ticks));
end

end
