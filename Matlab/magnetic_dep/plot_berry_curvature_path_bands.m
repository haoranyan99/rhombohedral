function H = plot_berry_curvature_path_bands(txt_file)


default_root = "E:\rg_master\rhombohedral\data";

band_list = [4, 5, 6];
energy_unit = "eV";
energy_shift_eV = 0.0;
color_energy_window_eV = [0.45, 0.9];
berry_clim = [-100,100];         % [] = auto; or set manually, e.g. [-10, 10]

FS = 15;
LW = 3.0;
save_figure = true;
out_png = "";

if nargin < 1 || strlength(string(txt_file)) == 0
    if ~isfolder(default_root)
        default_root = pwd;
    end

    [fn, fp] = uigetfile( ...
        fullfile(default_root, "berry_curvature_*.txt"), ...
        "Select berry_curvature path txt file");

    if isequal(fn, 0)
        H = struct();
        return;
    end

    txt_file = string(fullfile(fp, fn));
else
    txt_file = string(txt_file);
end

[data, columns, xtick_pos, xtick_lab, meta] = read_berry_txt_(txt_file);

if ~isfield(columns, "s")
    error("This script expects path output with column 's'. Selected file seems to be kmesh output.");
end

s = data(:, columns.s);
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
        eylabel = "E (eV)";
    case "mev"
        escale = 1e3;
        eylabel = "E (meV)";
    otherwise
        error('energy_unit must be "eV" or "meV".');
end

fig = figure("Color", "w", "Units", "pixels", "Position", [120 120 920 620]);
ax = axes(fig);
hold(ax, "on");
box(ax, "on");
set(ax, "FontSize", FS, "LineWidth", 1.3, "TickDir", "in", "Box", "on");
colormap(ax, turbo(256));

all_omega_colored = [];
all_omega_selected = [];
for b = bands
    ecol = columns.(sprintf("E_b%d", b));
    ocol = columns.(sprintf("berry_A2_b%d", b));
    E_eV = data(:, ecol) - energy_shift_eV;
    mask = E_eV >= color_energy_window_eV(1) ...
         & E_eV <= color_energy_window_eV(2);
    all_omega_colored = [all_omega_colored; data(mask, ocol)]; %#ok<AGROW>
    all_omega_selected = [all_omega_selected; data(:, ocol)]; %#ok<AGROW>
end
if isempty(all_omega_colored)
    warning("No k-path points fall inside color_energy_window_eV = [%g, %g]. Bands will be gray; colorbar uses all selected-band Berry values.", ...
        color_energy_window_eV(1), color_energy_window_eV(2));
    all_omega_colored = all_omega_selected;
end
if isempty(berry_clim)
    clim(ax, minmax_clim_(all_omega_colored, "minmax"));
else
    clim(ax, berry_clim);
end

H = struct();
H.file = txt_file;
H.bands = bands;
H.curves = struct([]);

for ib = 1:numel(bands)
    b = bands(ib);
    ecol = columns.(sprintf("E_b%d", b));
    ocol = columns.(sprintf("berry_A2_b%d", b));

    E_eV = data(:, ecol) - energy_shift_eV;
    E = E_eV * escale;
    omega = data(:, ocol);

    hline = plot(ax, s, E, ...
        "Color", [0.72 0.72 0.72], ...
        "LineWidth", LW, ...
        "DisplayName", sprintf("band %d", b));

    mask = E_eV >= color_energy_window_eV(1) ...
         & E_eV <= color_energy_window_eV(2);
    ranges = contiguous_ranges_(mask);

    for ir = 1:size(ranges, 1)
        idx = ranges(ir, 1):ranges(ir, 2);
        if numel(idx) < 2
            continue;
        end

        surface(ax, ...
            [s(idx) s(idx)], ...
            [E(idx) E(idx)], ...
            zeros(numel(idx), 2), ...
            [omega(idx) omega(idx)], ...
            "FaceColor", "none", ...
            "EdgeColor", "interp", ...
            "LineWidth", LW, ...
            "HandleVisibility", "off");
    end

    H.curves(ib).band = b;
    H.curves(ib).s = s;
    H.curves(ib).E_eV = E_eV;
    H.curves(ib).berry_A2 = omega;
    H.curves(ib).color_mask = mask;
    H.curves(ib).handle = hline;
end

if ~isempty(xtick_pos)
    xticks(ax, xtick_pos);
    xticklabels(ax, xtick_lab);
end

xlabel(ax, "k path");
ylabel(ax, eylabel);
cb = colorbar(ax);
cb.Label.String = "\Omega_z (Angstrom^2)";
cb.LineWidth = 1.0;
title(ax, make_title_(txt_file, meta), "Interpreter", "none");
legend(ax, "Location", "best", "Box", "off");
xlim(ax, [min(s), max(s)]);

if save_figure
    if strlength(out_png) == 0
        out_png = replace(txt_file, ".txt", "_berry_colored_bands.png");
    end
    exportgraphics(fig, out_png, "Resolution", 300);
    fprintf("Saved figure: %s\n", out_png);
end

H.figure = fig;
H.axes = ax;
H.colorbar = cb;
H.columns = columns;
H.xtick_pos = xtick_pos;
H.xtick_lab = xtick_lab;
H.meta = meta;
H.color_energy_window_eV = color_energy_window_eV;
H.berry_clim = clim(ax);

end

function [data, columns, xtick_pos, xtick_lab, meta] = read_berry_txt_(txt_file)

fid = fopen(txt_file, "r");
if fid < 0
    error("Cannot open file: %s", txt_file);
end
cleanup = onCleanup(@() fclose(fid));

columns_line = "";
xtick_pos = [];
xtick_lab = strings(0);
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
    elseif startsWith(body, "xtick =")
        tok = split(strtrim(extractAfter(body, "xtick =")));
        if numel(tok) >= 2
            xtick_pos(end+1, 1) = str2double(tok(1)); %#ok<AGROW>
            xtick_lab(end+1, 1) = strjoin(tok(2:end), " "); %#ok<AGROW>
        end
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

function clim_value = minmax_clim_(x, mode)

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
elseif mode == "symmetric"
    c = max(abs([xmin, xmax]));
    clim_value = [-c, c];
else
    clim_value = [xmin, xmax];
end

end

function ranges = contiguous_ranges_(mask)

mask = logical(mask(:));
d = diff([false; mask; false]);
ranges = [find(d == 1), find(d == -1) - 1];

end

function ttl = make_title_(txt_file, meta)

[~, base, ~] = fileparts(txt_file);
ttl = string(base);
if isfield(meta, "model")
    ttl = ttl + " | model=" + string(meta.model);
end
if isfield(meta, "Dfield_eV")
    ttl = ttl + " | D=" + string(meta.Dfield_eV) + " eV";
end
ttl = char(ttl);

end
