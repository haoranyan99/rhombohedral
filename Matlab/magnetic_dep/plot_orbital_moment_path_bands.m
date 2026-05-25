function H = plot_orbital_moment_path_bands(txt_file)

clc; close all;

% ============================================================
% USER SETTINGS
% ============================================================
default_root = "E:\rg_master\rhombohedral\data";

band_list = [4, 5, 6];      % [] means auto-detect all E_b*/m_orb_muB_b* pairs
energy_unit = "eV";      % "eV" or "meV"
energy_shift_eV = 0.0;   % subtract this before plotting
color_energy_window_eV = [0.5, 0.8];  % only this energy range is colored by m

FS = 15;
LW = 3.0;

save_figure = true;
out_png = "";            % "" means save next to the selected txt file

% ============================================================
% Select orbital moment path output
% ============================================================
if nargin < 1 || strlength(string(txt_file)) == 0
    if ~isfolder(default_root)
        default_root = pwd;
    end

    [fn, fp] = uigetfile( ...
        fullfile(default_root, "orbital_moment_*.txt"), ...
        "Select orbital_moment path txt file");

    if isequal(fn, 0)
        H = struct();
        return;
    end

    txt_file = string(fullfile(fp, fn));
else
    txt_file = string(txt_file);
end

[data, columns, xtick_pos, xtick_lab, meta] = read_orbital_moment_txt_(txt_file);

if ~isfield(columns, "s")
    error("This script expects path output with column 's'. Selected file seems to be kmesh output.");
end

s = data(:, columns.s);

if numel(unique(s)) < 2
    error("Path coordinate column 's' has too few distinct values.");
end

bands = detect_bands_(columns);
if isempty(bands)
    error("No E_b*/m_orb_muB_b* column pairs found.");
end

if ~isempty(band_list)
    bands = intersect(bands, double(band_list(:).'), "stable");
end

if isempty(bands)
    error("band_list does not match any E_b*/m_orb_muB_b* pair in this file.");
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

% ============================================================
% Figure
% ============================================================
fig = figure("Color", "w", "Units", "pixels", "Position", [120 120 920 620]);
ax = axes(fig);
hold(ax, "on");
box(ax, "on");

set(ax, ...
    "FontSize", FS, ...
    "LineWidth", 1.3, ...
    "TickDir", "in", ...
    "Box", "on");

cmap = turbo(256);
colormap(ax, cmap);

all_m_colored = [];
for b = bands
    ecol = columns.(sprintf("E_b%d", b));
    mcol = columns.(sprintf("m_orb_muB_b%d", b));
    E_eV = data(:, ecol) - energy_shift_eV;
    mask = E_eV >= color_energy_window_eV(1) ...
         & E_eV <= color_energy_window_eV(2);
    all_m_colored = [all_m_colored; data(mask, mcol)]; %#ok<AGROW>
end

mmin = min(all_m_colored, [], "omitnan");
mmax = max(all_m_colored, [], "omitnan");
if ~(isfinite(mmin) && isfinite(mmax))
    mmin = -1.0;
    mmax = 1.0;
elseif mmin == mmax
    delta = max(abs(mmin), 1.0) * 1e-6;
    mmin = mmin - delta;
    mmax = mmax + delta;
end
clim(ax, [mmin, mmax]);

H = struct();
H.file = txt_file;
H.bands = bands;
H.curves = struct([]);

for ib = 1:numel(bands)
    b = bands(ib);

    ecol = columns.(sprintf("E_b%d", b));
    mcol = columns.(sprintf("m_orb_muB_b%d", b));

    E_eV = data(:, ecol) - energy_shift_eV;
    E = E_eV * escale;
    m = data(:, mcol);

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

        handle_visibility = "off";
        if ir == 1
            handle_visibility = "off";
        end

        surface(ax, ...
            [s(idx) s(idx)], ...
            [E(idx) E(idx)], ...
            zeros(numel(idx), 2), ...
            [m(idx) m(idx)], ...
            "FaceColor", "none", ...
            "EdgeColor", "interp", ...
            "LineWidth", LW, ...
            "HandleVisibility", handle_visibility);
    end

    H.curves(ib).band = b;
    H.curves(ib).s = s;
    H.curves(ib).E_eV = E_eV;
    H.curves(ib).E_plot = E;
    H.curves(ib).m_orb_muB = m;
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
cb.Label.String = "m_{orb} (\mu_B)";
cb.LineWidth = 1.0;

title(ax, make_title_(txt_file, meta), "Interpreter", "none");
legend(ax, "Location", "best", "Box", "off");

xlim(ax, [min(s), max(s)]);

if save_figure
    if strlength(out_png) == 0
        out_png = replace(txt_file, ".txt", "_colored_bands.png");
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

end

function [data, columns, xtick_pos, xtick_lab, meta] = read_orbital_moment_txt_(txt_file)

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
    mfield = sprintf("m_orb_muB_b%d", b);
    if isfield(columns, mfield)
        bands(end+1) = b; %#ok<AGROW>
    end
end

bands = sort(unique(bands));

end

function ttl = make_title_(txt_file, meta)

[~, base, ~] = fileparts(txt_file);
parts = string(base);

if isfield(meta, "model")
    parts = parts + " | model=" + string(meta.model);
end

if isfield(meta, "Dfield_eV")
    parts = parts + " | D=" + string(meta.Dfield_eV) + " eV";
end

ttl = char(parts);

end

function ranges = contiguous_ranges_(mask)

mask = logical(mask(:));
if isempty(mask)
    ranges = zeros(0, 2);
    return;
end

d = diff([false; mask; false]);
starts = find(d == 1);
stops = find(d == -1) - 1;
ranges = [starts, stops];

end
