function H = plot_orbital_moment_kmesh_maps(txt_file)

% ============================================================
% USER SETTINGS
% ============================================================
default_root = "E:\rg_master\rhombohedral\data";

band_list = [4, 5];      % [] means auto-detect all E_b*/m_orb_muB_b* pairs
energy_unit = "eV";      % "eV" or "meV"
energy_shift_eV = 0.0;   % subtract this before plotting

xy_mode = "kxy";         % "kxy" or "ij"
marker_size = 12;

m_clim_mode = "minmax";  % "minmax" or "symmetric"

FS = 15;

save_figure = true;
out_dir = "";            % "" means save next to the selected txt file

% ============================================================
% Select orbital moment kmesh output
% ============================================================
if nargin < 1 || strlength(string(txt_file)) == 0
    if ~isfolder(default_root)
        default_root = pwd;
    end

    [fn, fp] = uigetfile( ...
        fullfile(default_root, "orbital_moment_*.txt"), ...
        "Select orbital_moment kmesh txt file");

    if isequal(fn, 0)
        H = struct();
        return;
    end

    txt_file = string(fullfile(fp, fn));
else
    txt_file = string(txt_file);
end

[data, columns, meta] = read_orbital_moment_txt_(txt_file);

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

for ib = 1:numel(bands)
    b = bands(ib);

    ecol = columns.(sprintf("E_b%d", b));
    mcol = columns.(sprintf("m_orb_muB_b%d", b));

    E = (data(:, ecol) - energy_shift_eV) * escale;
    m = data(:, mcol);

    figE = make_map_figure_( ...
        x, y, E, ...
        sprintf("band %d: E(k)", b), ...
        e_label, ...
        x_label, y_label, ...
        parula(256), ...
        [], ...
        FS, marker_size);

    m_clim = choose_m_clim_(m, m_clim_mode);
    figM = make_map_figure_( ...
        x, y, m, ...
        sprintf("band %d: m_{orb}(k)", b), ...
        "m_{orb} (\mu_B)", ...
        x_label, y_label, ...
        turbo(256), ...
        m_clim, ...
        FS, marker_size);

    if save_figure
        outE = fullfile(out_dir, sprintf("%s_band%d_E_%s_map.png", base, b, e_tag));
        outM = fullfile(out_dir, sprintf("%s_band%d_morb_muB_map.png", base, b));

        exportgraphics(figE, outE, "Resolution", 300);
        exportgraphics(figM, outM, "Resolution", 300);

        fprintf("Saved figure: %s\n", outE);
        fprintf("Saved figure: %s\n", outM);
    end

    H.maps(ib).band = b;
    H.maps(ib).E_plot = E;
    H.maps(ib).m_orb_muB = m;
    H.maps(ib).figure_E = figE;
    H.maps(ib).figure_m = figM;
end

end

function fig = make_map_figure_(x, y, val, ttl, cb_label, x_label, y_label, cmap, clim_value, FS, marker_size)

fig = figure("Color", "w", "Units", "pixels", "Position", [120 120 720 620]);
ax = axes(fig);

scatter(ax, x, y, marker_size, val, "filled", ...
    "MarkerEdgeAlpha", 0.0, ...
    "MarkerFaceAlpha", 1.0);

axis(ax, "equal");
box(ax, "on");
set(ax, ...
    "FontSize", FS, ...
    "LineWidth", 1.3, ...
    "TickDir", "in", ...
    "Box", "on");

colormap(ax, cmap);

if ~isempty(clim_value)
    clim(ax, clim_value);
end

cb = colorbar(ax);
cb.Label.String = cb_label;
cb.LineWidth = 1.0;

xlabel(ax, x_label);
ylabel(ax, y_label);
title(ax, ttl, "Interpreter", "none");

end

function clim_value = choose_m_clim_(m, mode)

mmin = min(m, [], "omitnan");
mmax = max(m, [], "omitnan");

if ~(isfinite(mmin) && isfinite(mmax))
    clim_value = [-1, 1];
    return;
end

if mmin == mmax
    delta = max(abs(mmin), 1.0) * 1e-6;
    clim_value = [mmin - delta, mmax + delta];
    return;
end

switch lower(mode)
    case "minmax"
        clim_value = [mmin, mmax];
    case "symmetric"
        c = max(abs([mmin, mmax]));
        clim_value = [-c, c];
    otherwise
        error('m_clim_mode must be "minmax" or "symmetric".');
end

end

function [data, columns, meta] = read_orbital_moment_txt_(txt_file)

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
    mfield = sprintf("m_orb_muB_b%d", b);
    if isfield(columns, mfield)
        bands(end+1) = b; %#ok<AGROW>
    end
end

bands = sort(unique(bands));

end
