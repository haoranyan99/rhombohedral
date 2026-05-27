function H = plot_berry_curvature_kmesh_maps(txt_file)

data_folder = "E:\rg_master\rhombohedral\data\berry_curvature_sk_50_b1b2_0.00125\D-0.010";

band_list = 4;
energy_unit = "eV";
energy_shift_eV = 0.0;
xy_mode = "kxy";          % "kxy" or "ij"
omega_clim_mode = "minmax"; % "minmax" or "symmetric"
berry_clim = [];          % [] = auto; or set manually, e.g. [-10, 10]
marker_size = 12;
FS = 15;
save_figure = true;
save_fig_file = false;
out_dir = "";

if nargin < 1 || strlength(string(txt_file)) == 0
    txt_file = latest_txt_in_folder_(data_folder, "berry_curvature_*.txt");

    if strlength(txt_file) == 0
        error("Cannot find berry_curvature_*.txt under folder: %s", data_folder);
    end

    fprintf("Using file: %s\n", txt_file);
else
    txt_file = string(txt_file);
end

[data, columns, meta] = read_berry_txt_(txt_file);

if isfield(columns, "s")
    error("This script expects 2D kmesh output. Selected file seems to be kpath output.");
end

required = ["kx", "ky", "iq", "jq"];
for r = required
    if ~isfield(columns, r)
        error("Missing required column '%s'.", r);
    end
end

bands = detect_bands_(columns);
if ~isempty(band_list)
    bands = intersect(bands, double(band_list(:).'), "stable");
end
if isempty(bands)
    error("No requested E_b*/berry_A2_b* column pairs found.");
end

switch lower(energy_unit)
    case "ev"
        escale = 1.0;
        e_label = "E (eV)";
        e_tag = "eV";
    case "mev"
        escale = 1e3;
        e_label = "E (meV)";
        e_tag = "meV";
    otherwise
        error('energy_unit must be "eV" or "meV".');
end

switch lower(xy_mode)
    case "kxy"
        x = data(:, columns.kx);
        y = data(:, columns.ky);
        x_label = "k_x (Angstrom^{-1})";
        y_label = "k_y (Angstrom^{-1})";
    case "ij"
        x = data(:, columns.iq);
        y = data(:, columns.jq);
        x_label = "i_q";
        y_label = "j_q";
    otherwise
        error('xy_mode must be "kxy" or "ij".');
end

if strlength(out_dir) == 0
    out_dir = string(fileparts(txt_file));
else
    out_dir = string(out_dir);
end
if save_figure && ~isfolder(out_dir)
    mkdir(out_dir);
end

[~, base, ~] = fileparts(txt_file);

H = struct();
H.file = txt_file;
H.bands = bands;
H.maps = struct([]);
H.meta = meta;
H.columns = columns;
H.berry_clim = berry_clim;

for ib = 1:numel(bands)
    b = bands(ib);
    ecol = columns.(sprintf("E_b%d", b));
    ocol = columns.(sprintf("berry_A2_b%d", b));

    E = (data(:, ecol) - energy_shift_eV) * escale;
    omega = data(:, ocol);

    figPair = make_pair_map_figure_( ...
        x, y, E, ...
        omega, ...
        sprintf("band %d", b), ...
        e_label, "\Omega_z (Angstrom^2)", ...
        x_label, y_label, ...
        parula(256), turbo(256), ...
        [], ...
        choose_clim_(omega, omega_clim_mode, berry_clim), ...
        FS, marker_size);

    if save_figure
        outPair = fullfile(out_dir, sprintf("%s_band%d_E_%s_berry_A2_map.png", base, b, e_tag));
        exportgraphics(figPair, outPair, "Resolution", 300);
        if save_fig_file
            savefig(figPair, replace(outPair, ".png", ".fig"));
        end
        fprintf("Saved figure: %s\n", outPair);
    end

    H.maps(ib).band = b;
    H.maps(ib).E_plot = E;
    H.maps(ib).berry_A2 = omega;
    H.maps(ib).figure_pair = figPair;
end

end

function f = latest_txt_in_folder_(folder, pattern)
d = dir(fullfile(folder, pattern));
d = d(~[d.isdir]);
if isempty(d)
    f = "";
    return;
end
[~, idx] = max([d.datenum]);
f = string(fullfile(d(idx).folder, d(idx).name));
end

function fig = make_pair_map_figure_(x, y, val_left, val_right, ttl, ...
    cb_label_left, cb_label_right, x_label, y_label, cmap_left, cmap_right, ...
    clim_left, clim_right, FS, marker_size)

fig = figure("Color", "w", "Units", "pixels", "Position", [120 120 1160 540]);
tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");

ax1 = nexttile;
scatter(ax1, x, y, marker_size, val_left, "filled", ...
    "MarkerEdgeAlpha", 0.0, ...
    "MarkerFaceAlpha", 1.0);
axis(ax1, "equal");
box(ax1, "on");
set(ax1, "FontSize", FS, "LineWidth", 1.3, "TickDir", "in", "Box", "on");
colormap(ax1, cmap_left);
if ~isempty(clim_left)
    clim(ax1, clim_left);
end
cb1 = colorbar(ax1);
cb1.Label.String = cb_label_left;
cb1.LineWidth = 1.0;
xlabel(ax1, x_label);
ylabel(ax1, y_label);
title(ax1, "energy", "Interpreter", "none");

ax2 = nexttile;
scatter(ax2, x, y, marker_size, val_right, "filled", ...
    "MarkerEdgeAlpha", 0.0, ...
    "MarkerFaceAlpha", 1.0);
axis(ax2, "equal");
box(ax2, "on");
set(ax2, "FontSize", FS, "LineWidth", 1.3, "TickDir", "in", "Box", "on");
colormap(ax2, cmap_right);
if ~isempty(clim_right)
    clim(ax2, clim_right);
end
cb2 = colorbar(ax2);
cb2.Label.String = cb_label_right;
cb2.LineWidth = 1.0;
xlabel(ax2, x_label);
ylabel(ax2, y_label);
title(ax2, "Berry curvature", "Interpreter", "none");

sgtitle(fig, ttl, "Interpreter", "none", "FontSize", FS);

end

function clim_value = choose_clim_(x, mode, manual_clim)

if ~isempty(manual_clim)
    clim_value = manual_clim;
    return;
end

x = x(isfinite(x));
if isempty(x)
    clim_value = [-1, 1];
    return;
end

xmin = min(x, [], "omitnan");
xmax = max(x, [], "omitnan");
if ~(isfinite(xmin) && isfinite(xmax))
    clim_value = [-1, 1];
elseif xmin == xmax
    d = max(abs(xmin), 1.0) * 1e-6;
    clim_value = [xmin - d, xmax + d];
elseif lower(mode) == "symmetric"
    c = max(abs([xmin, xmax]));
    clim_value = [-c, c];
else
    clim_value = [xmin, xmax];
end

end

function [data, columns, meta] = read_berry_txt_(txt_file)

fid = fopen(txt_file, "r");
if fid < 0
    error("Cannot open file: %s", txt_file);
end
cleanup = onCleanup(@() fclose(fid));

columns_line = "";
meta = struct();
while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    line = string(strtrim(line));
    if ~startsWith(line, "#")
        continue;
    end
    body = strtrim(extractAfter(line, 1));
    if startsWith(body, "columns =")
        columns_line = strtrim(extractAfter(body, "columns ="));
    elseif contains(body, "=")
        ieq = strfind(body, "=");
        ieq = ieq(1);
        key = matlab.lang.makeValidName(strtrim(extractBefore(body, ieq)));
        val = strtrim(extractAfter(body, ieq));
        meta.(key) = val;
    end
end

if strlength(columns_line) == 0
    error("Cannot find '# columns =' header in %s", txt_file);
end

names = split(columns_line);
columns = struct();
for i = 1:numel(names)
    columns.(matlab.lang.makeValidName(names(i))) = i;
end

data = readmatrix(txt_file, "FileType", "text", "CommentStyle", "#");
data = data(all(isfinite(data), 2), :);
if isempty(data)
    error("No numeric data read from %s", txt_file);
end

end

function bands = detect_bands_(columns)

fields = string(fieldnames(columns));
bands = [];
for f = fields(:).'
    tok = regexp(f, "^E_b(\d+)$", "tokens", "once");
    if isempty(tok)
        continue;
    end
    b = str2double(tok{1});
    if isfield(columns, sprintf("berry_A2_b%d", b))
        bands(end+1) = b; %#ok<AGROW>
    end
end
bands = sort(unique(bands));

end
