function OUT = plot_hyst_phasemap_polar_doping_from_mat()

clc; close all;

% =========================
% USER SETTINGS
% =========================
T_target = 6.4;          % fixed temperature, K
T_tol = 1e-8;

quantity = "psi";        % "psi" or "X"

FS = 16;
doping_N = 300;

use_abs_diff = true;

% Manual display range. Use [] to keep the full common valid range.
xlim_plot = [-2.5, 0.1];   % doping (10^12 cm^-2)
ylim_plot = [0, 2.5];    % artificial polar (meV)

hot_yellow_value_frac = 0.5;
hot_yellow_color_pos = 0.5;
hot_white_value = 1.2;
hot_white_color_pos = 0.06;
hot_red_value = 4.1;
hot_red_color_pos = 0.65;
white_pos = 0.75;        % used only for signed fwd-bwd maps

default_dir = fullfile(pwd, "Matlab", "hyst_data");

% =========================
% Load MAT
% =========================
if ~isfolder(default_dir)
    default_dir = pwd;
end

[fname, fpath] = uigetfile(fullfile(default_dir, "*.mat"), ...
    "Select hyst_allT_allArtificialPolar.mat");

if isequal(fname, 0)
    OUT = [];
    return;
end

matfile = fullfile(fpath, fname);
S = load(matfile);

if ~isfield(S, "HYST")
    error("Selected MAT file does not contain variable HYST: %s", matfile);
end

HYST = S.HYST(:);

valid = arrayfun(@(x) scalar_field_(x, "T") && ...
    scalar_field_(x, "artificial_polar"), HYST);
entries = HYST(valid);

if isempty(entries)
    error("No valid HYST entries found in %s.", matfile);
end

% =========================
% Select fixed temperature
% =========================
E = entries(abs([entries.T] - T_target) <= T_tol);

if isempty(E)
    Tall = unique([entries.T]);
    error("No HYST entries found for T=%.6g K. Available T: %s", ...
        T_target, mat2str(Tall));
end

% If the same polar was saved more than once, keep the last one in HYST order.
[~, keep_idx] = unique([E.artificial_polar], "last");
E = E(keep_idx);

% =========================
% Build common doping grid
% =========================
dop_min = -inf;
dop_max = inf;

for i = 1:numel(E)
    [xf, ~, xb, ~] = get_quantity_curves_(E(i), quantity);

    xf = xf(isfinite(xf));
    xb = xb(isfinite(xb));

    if isempty(xf) || isempty(xb)
        continue;
    end

    dop_min = max(dop_min, max(min(xf), min(xb)));
    dop_max = min(dop_max, min(max(xf), max(xb)));
end

if ~isempty(xlim_plot)
    dop_min = max(dop_min, min(xlim_plot));
    dop_max = min(dop_max, max(xlim_plot));
end

if ~(isfinite(dop_min) && isfinite(dop_max) && dop_max > dop_min)
    error("Cannot build a common doping range for T=%.6g K.", T_target);
end

dop_grid = linspace(dop_min, dop_max, doping_N);

% =========================
% Compute Z(polar,doping)
% =========================
polar = nan(numel(E), 1);
Z = nan(numel(E), doping_N);

for i = 1:numel(E)
    R = E(i);
    polar(i) = R.artificial_polar;

    [xf, yf, xb, yb] = get_quantity_curves_(R, quantity);

    okf = isfinite(xf) & isfinite(yf);
    okb = isfinite(xb) & isfinite(yb);

    xf = xf(okf);
    yf = yf(okf);
    xb = xb(okb);
    yb = yb(okb);

    if numel(xf) < 3 || numel(xb) < 3
        continue;
    end

    [xf, ord] = sort(xf); yf = yf(ord);
    [xb, ord] = sort(xb); yb = yb(ord);

    [xf, ia] = unique(xf); yf = yf(ia);
    [xb, ia] = unique(xb); yb = yb(ia);

    yf_grid = interp1(xf, yf, dop_grid, "linear", NaN);
    yb_grid = interp1(xb, yb, dop_grid, "linear", NaN);

    if use_abs_diff
        Z(i,:) = abs(yf_grid - yb_grid);
    else
        Z(i,:) = yf_grid - yb_grid;
    end
end

okrow = isfinite(polar) & any(isfinite(Z), 2);
polar = polar(okrow);
Z = Z(okrow,:);

[polar, ord] = sort(polar);
Z = Z(ord,:);

if isempty(polar)
    error("No valid polar-doping hysteresis map data found.");
end

% =========================
% Restrict display range
% =========================
maskP = true(size(polar));
if ~isempty(ylim_plot)
    maskP = polar >= min(ylim_plot) & polar <= max(ylim_plot);
end

polar_show = polar(maskP);
Z_show = Z(maskP,:);
dop_show = dop_grid;

if isempty(polar_show) || isempty(dop_show) || all(~isfinite(Z_show(:)))
    error("No valid data inside the selected plot range.");
end

% =========================
% Color range from displayed data
% =========================
if use_abs_diff
    zmin = 0;
else
    zmin = min(Z_show(:), [], "omitnan");
end

zmax = max(Z_show(:), [], "omitnan");

if ~(zmax > zmin)
    zmax = zmin + eps;
end

% =========================
% Plot
% =========================
fig = figure("Color","w","Position",[120 120 850 620]);
ax = axes(fig);

imagesc(ax, dop_show, polar_show, Z_show);
set(ax, "YDir", "normal");

if use_abs_diff
    zred = min(max(hot_red_value, zmin), zmax);
    zyellow_auto = zmin + hot_yellow_value_frac * (zmax - zmin);
    zyellow = min(zyellow_auto, hot_white_value + 0.5 * (zred - hot_white_value));
    zyellow = min(max(zyellow, zmin), zmax);
    white_value_frac = (hot_white_value - zmin) / (zmax - zmin);
    yellow_value_frac = (zyellow - zmin) / (zmax - zmin);
    red_value_frac = (zred - zmin) / (zmax - zmin);
    positive_colormap = positive_hot_colormap_( ...
        256, white_value_frac, hot_white_color_pos, ...
        yellow_value_frac, hot_yellow_color_pos, ...
        red_value_frac, hot_red_color_pos);
    clim(ax, [zmin, zmax]);
else
    zyellow = NaN;
    clim(ax, [zmin, zmax]);
end

cb = colorbar(ax);
cb.TickLabelInterpreter = "latex";
cb.TickDirection = "in";
cb.LineWidth = 1.2;

if use_abs_diff
    colormap(ax, positive_colormap);
    [tick_pos, tick_label] = positive_hot_colorbar_ticks_( ...
        zmin, zyellow, zmax);
    cb.Ticks = tick_pos;
    cb.TickLabels = tick_label;
    cb.Label.String = sprintf("$||%s|_{\\rm fwd}-|%s|_{\\rm bwd}|$", ...
        quantity, quantity);
else
    colormap(ax, red_white_blue_skewed_(256, white_pos));
    cb.Label.String = sprintf("$|%s|_{\\rm fwd}-|%s|_{\\rm bwd}$", ...
        quantity, quantity);
end
cb.Label.Interpreter = "latex";

xlabel(ax, "Doping $(10^{12}\,\mathrm{cm}^{-2})$", "Interpreter","latex");
ylabel(ax, "Artificial polar / magnetic field proxy (meV)", "Interpreter","latex");

title(ax, sprintf("$T=%.3f\\,\\mathrm{K}$, polar-doping hysteresis map of $|%s|$", ...
    T_target, quantity), "Interpreter","latex");

set(ax, ...
    "FontSize", FS, ...
    "TickLabelInterpreter","latex", ...
    "LineWidth", 1.3, ...
    "TickDir","in", ...
    "Box","on");

if ~isempty(xlim_plot)
    xlim(ax, [min(xlim_plot), max(xlim_plot)]);
end
if ~isempty(ylim_plot)
    ylim(ax, [min(ylim_plot), max(ylim_plot)]);
end
grid(ax, "off");

OUT = struct();
OUT.matfile = string(matfile);
OUT.T_target = T_target;
OUT.T_tol = T_tol;
OUT.quantity = quantity;
OUT.use_abs_diff = use_abs_diff;
OUT.doping = dop_show(:);
OUT.polar = polar_show(:);
OUT.Z = Z_show;
OUT.color_range = [zmin, zmax];
OUT.fig = fig;

fprintf("\nLoaded: %s\n", matfile);
fprintf("T = %.6g K | quantity = %s\n", T_target, quantity);
fprintf("Displayed doping range = [%.6g, %.6g]\n", dop_show(1), dop_show(end));
fprintf("Displayed polar range = [%.6g, %.6g]\n", min(polar_show), max(polar_show));
fprintf("Color range from displayed data = [%.6g, %.6g]\n", zmin, zmax);

end

function [xf, yf, xb, yb] = get_quantity_curves_(R, quantity)

xf = R.dop_fwd(:);
xb = R.dop_bwd(:);

switch lower(string(quantity))
    case "psi"
        if isfield(R, "psi_f_plot") && isfield(R, "psi_b_plot")
            yf = R.psi_f_plot(:);
            yb = R.psi_b_plot(:);
        elseif isfield(R, "abspsi_fwd") && isfield(R, "abspsi_bwd")
            yf = R.abspsi_fwd(:);
            yb = R.abspsi_bwd(:);
        elseif isfield(R, "psi_fwd") && isfield(R, "psi_bwd")
            yf = abs(R.psi_fwd(:));
            yb = abs(R.psi_bwd(:));
        else
            error("Cannot find psi data in HYST entry.");
        end

    case "x"
        if isfield(R, "X_f_plot") && isfield(R, "X_b_plot")
            yf = R.X_f_plot(:);
            yb = R.X_b_plot(:);
        elseif isfield(R, "absX_fwd") && isfield(R, "absX_bwd")
            yf = R.absX_fwd(:);
            yb = R.absX_bwd(:);
        elseif isfield(R, "X_fwd") && isfield(R, "X_bwd")
            yf = abs(R.X_fwd(:));
            yb = abs(R.X_bwd(:));
        else
            error("Cannot find X data in HYST entry.");
        end

    otherwise
        error("quantity must be psi or X.");
end

end

function tf = scalar_field_(R, fieldname)

tf = isfield(R, fieldname) && isscalar(R.(fieldname)) && ...
    isnumeric(R.(fieldname)) && isfinite(R.(fieldname));

end

function cmap = red_white_blue_skewed_(N, white_pos)

if nargin < 1
    N = 256;
end
if nargin < 2
    white_pos = 0.5;
end

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

function cmap = positive_hot_colormap_(N, white_value_frac, white_color_pos, yellow_value_frac, yellow_color_pos, red_value_frac, red_color_pos)

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

function [ticks, labels] = positive_hot_colorbar_ticks_(zmin, zyellow, zmax)

ticks = unique([ceil(zmin), round(zyellow), floor(zmax)]);
ticks = ticks(ticks >= zmin & ticks <= zmax);

if numel(ticks) < 3
    ticks = unique([zmin, zyellow, zmax]);
    labels = cellstr(compose("%.4g", ticks));
else
    labels = cellstr(compose("%.0f", ticks));
end

end
