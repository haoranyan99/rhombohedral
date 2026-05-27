function H = plot_orbital_moment_kmesh_maps(txt_file)

% ============================================================
% USER SETTINGS
% ============================================================
data_folder = "E:\rg_master\rhombohedral\data\orbital_moment_sk_50_b1b2_0.00125\D-0.010";

band_list = 4;         % choose one band, or a short list; [] means all bands
energy_unit = "eV";      % "eV" or "meV"
energy_shift_eV = 0.0;   % subtract this before plotting
B_list_T = [1, 5];       % plot orbital energy shifts for these fields
valley_choice = "plus";  % "plus" or "minus"
spin_choice = "up";      % "up" or "down"
g_factor = 2.0;
plot_energy_map = false;

xy_mode = "kxy";         % "kxy" or "ij"
marker_size = 12;

m_clim_mode = "minmax";  % "minmax" or "symmetric"

FS = 15;

save_figure = true;
save_fig_file = false;
out_dir = "";            % "" means save next to the selected txt file

% ============================================================
% Select orbital moment kmesh output
% ============================================================
if nargin < 1 || strlength(string(txt_file)) == 0
    txt_file = latest_txt_in_folder_(data_folder, "orbital_moment_*.txt");

    if strlength(txt_file) == 0
        error("Cannot find orbital_moment_*.txt under folder: %s", data_folder);
    end

    fprintf("Using file: %s\n", txt_file);
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
H.B_list_T = B_list_T;
H.valley_choice = valley_choice;
H.spin_choice = spin_choice;
H.g_factor = g_factor;

tau = valley_sign_(valley_choice);
s_spin = spin_sign_(spin_choice);

for ib = 1:numel(bands)
    b = bands(ib);

    ecol = columns.(sprintf("E_b%d", b));
    mcol = columns.(sprintf("m_orb_muB_b%d", b));

    E = (data(:, ecol) - energy_shift_eV) * escale;
    m = data(:, mcol);
    m_valley = tau * m;

    if plot_energy_map
        figE = make_map_figure_( ...
            x, y, E, ...
            sprintf("band %d: E(k)", b), ...
            e_label, ...
            x_label, y_label, ...
            parula(256), ...
            [], ...
            FS, marker_size);

        if save_figure
            outE = fullfile(out_dir, sprintf("%s_band%d_E_%s_map.png", base, b, e_tag));
            exportgraphics(figE, outE, "Resolution", 300);
            fprintf("Saved figure: %s\n", outE);
        end
    else
        figE = [];
    end

    for iB = 1:numel(B_list_T)
        B_T = B_list_T(iB);
        dE_orb_meV = -m_valley * 5.7883818060e-2 * B_T;
        dE_total_meV = dE_orb_meV ...
            + s_spin * g_factor * 5.7883818060e-2 * B_T;

        m_clim = choose_m_clim_(m_valley, m_clim_mode);
        dE_clim = choose_m_clim_(dE_total_meV, "minmax");

        figPair = make_pair_map_figure_( ...
            x, y, ...
            m_valley, dE_total_meV, ...
            sprintf("band %d, valley %s, spin %s, B=%g T", ...
                b, valley_choice, spin_choice, B_T), ...
            "m_{orb} (\mu_B)", "\DeltaE_{orb+Z} (meV)", ...
            x_label, y_label, ...
            turbo(256), parula(256), ...
            m_clim, dE_clim, ...
            FS, marker_size);

        if save_figure
            btag = B_tag_(B_T);
            outPair = fullfile(out_dir, sprintf( ...
                "%s_band%d_valley_%s_spin_%s_morb_dE_%s_map.png", ...
                base, b, valley_choice, spin_choice, btag));
            exportgraphics(figPair, outPair, "Resolution", 300);
            if save_fig_file
                savefig(figPair, replace(outPair, ".png", ".fig"));
            end
            fprintf("Saved figure: %s\n", outPair);
        end

        H.maps(ib).B(iB).B_T = B_T;
        H.maps(ib).B(iB).m_orb_valley_muB = m_valley;
        H.maps(ib).B(iB).dE_orb_meV = dE_orb_meV;
        H.maps(ib).B(iB).dE_total_meV = dE_total_meV;
        H.maps(ib).B(iB).figure_pair = figPair;
    end

    H.maps(ib).band = b;
    H.maps(ib).E_plot = E;
    H.maps(ib).m_orb_muB = m_valley;
    H.maps(ib).figure_E = figE;
end

function tag = B_tag_(B_T)
s = sprintf("%.10f", B_T);
s = regexprep(s, "0+$", "");
s = regexprep(s, "\.$", "");
if s == "-0"
    s = "0";
end
tag = "B" + string(s) + "T";
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

end

function tau = valley_sign_(valley_choice)
switch lower(string(valley_choice))
    case {"plus", "+", "valley_plus"}
        tau = +1;
    case {"minus", "-", "valley_minus"}
        tau = -1;
    otherwise
        error('valley_choice must be "plus" or "minus".');
end
end

function s = spin_sign_(spin_choice)
switch lower(string(spin_choice))
    case {"up", "+", "spin_up"}
        s = +1;
    case {"down", "dn", "dwn", "-", "spin_down"}
        s = -1;
    otherwise
        error('spin_choice must be "up" or "down".');
end
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
title(ax1, "orbital moment", "Interpreter", "none");

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
title(ax2, "energy shift", "Interpreter", "none");

sgtitle(fig, ttl, "Interpreter", "none", "FontSize", FS);

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
