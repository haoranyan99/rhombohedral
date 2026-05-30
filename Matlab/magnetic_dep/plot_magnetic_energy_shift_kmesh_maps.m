function H = plot_magnetic_energy_shift_kmesh_maps(txt_file)

% ============================================================
% USER SETTINGS
% ============================================================
data_folder = "E:\rg_master\rhombohedral\data\orbital_moment_sk_50_b1b2_0.00125\D-0.010";

band_list = 4;              % choose one band, a list, or [] for all
B_list_T = [1];
valley_choice = "minus";     % "plus", "minus", or "average"
spin_choice = "average";    % "up", "down", or "average"
valley_choice = "average"; 
spin_choice = "up";
g_factor = 2.0;

xy_mode = "kxy";            % "kxy" or "ij"
marker_size = 12;
shift_unit = "meV";         % "meV" or "eV"
shift_clim = [];            % [] = auto, or e.g. [-0.5, 0.5]
shift_clim_mode = "minmax"; % "minmax" or "symmetric"

FS = 15;
save_figure = true;
save_fig_file = false;
out_dir = "";              % "" means save next to selected txt

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
if ~isempty(band_list)
    bands = intersect(bands, double(band_list(:).'), "stable");
end
if isempty(bands)
    error("No requested E_b*/m_orb_muB_b* column pairs found.");
end

switch lower(shift_unit)
    case "mev"
        unit_scale = 1.0;
        cb_label = "\Delta E (meV)";
        unit_tag = "meV";
    case "ev"
        unit_scale = 1e-3;
        cb_label = "\Delta E (eV)";
        unit_tag = "eV";
    otherwise
        error('shift_unit must be "meV" or "eV".');
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
H.meta = meta;
H.columns = columns;
H.B_list_T = B_list_T;
H.valley_choice = valley_choice;
H.spin_choice = spin_choice;
H.g_factor = g_factor;
H.maps = struct([]);

for ib = 1:numel(bands)
    b = bands(ib);
    mcol = columns.(sprintf("m_orb_muB_b%d", b));
    m_plus_valley_muB = data(:, mcol);

    H.maps(ib).band = b;
    H.maps(ib).B = struct([]);

    for iB = 1:numel(B_list_T)
        B_T = B_list_T(iB);
        dE_meV = branch_average_shift_meV_( ...
            m_plus_valley_muB, B_T, g_factor, valley_choice, spin_choice);
        dE_plot = dE_meV * unit_scale;

        clim_value = choose_clim_(dE_plot, shift_clim_mode, shift_clim);

        ttl = sprintf("band %d, B=%g T, valley=%s, spin=%s", ...
            b, B_T, valley_choice, spin_choice);

        fig = make_shift_map_figure_( ...
            x, y, dE_plot, ttl, cb_label, ...
            x_label, y_label, turbo(256), clim_value, FS, marker_size);

        if save_figure
            btag = B_tag_(B_T);
            out_png = fullfile(out_dir, sprintf( ...
                "%s_band%d_energy_shift_valley_%s_spin_%s_%s_%s_map.png", ...
                base, b, valley_choice, spin_choice, btag, unit_tag));
            exportgraphics(fig, out_png, "Resolution", 300);
            if save_fig_file
                savefig(fig, replace(out_png, ".png", ".fig"));
            end
            fprintf("Saved figure: %s\n", out_png);
        end

        H.maps(ib).B(iB).B_T = B_T;
        H.maps(ib).B(iB).dE_meV = dE_meV;
        H.maps(ib).B(iB).dE_plot = dE_plot;
        H.maps(ib).B(iB).figure = fig;
    end
end

end

function dE_avg_meV = branch_average_shift_meV_(m_plus_valley_muB, B_T, g_factor, valley_choice, spin_choice)

muB_meV_per_T = 5.7883818060e-2;
taus = valley_choices_(valley_choice);
spins = spin_choices_(spin_choice);

dE_sum = zeros(size(m_plus_valley_muB));
n = 0;

for it = 1:numel(taus)
    for is = 1:numel(spins)
        tau = taus(it);
        s = spins(is);

        dE_orb = -tau * m_plus_valley_muB * muB_meV_per_T * B_T;
        dE_z = s * g_factor * muB_meV_per_T * B_T;

        dE_sum = dE_sum + dE_orb + dE_z;
        n = n + 1;
    end
end

dE_avg_meV = dE_sum / n;

end

function taus = valley_choices_(valley_choice)

switch lower(string(valley_choice))
    case {"plus", "+", "valley_plus"}
        taus = +1;
    case {"minus", "-", "valley_minus"}
        taus = -1;
    case {"average", "avg", "mean"}
        taus = [+1, -1];
    otherwise
        error('valley_choice must be "plus", "minus", or "average".');
end

end

function spins = spin_choices_(spin_choice)

switch lower(string(spin_choice))
    case {"up", "+", "spin_up"}
        spins = +1;
    case {"down", "dn", "dwn", "-", "spin_down"}
        spins = -1;
    case {"average", "avg", "mean"}
        spins = [+1, -1];
    otherwise
        error('spin_choice must be "up", "down", or "average".');
end

end

function fig = make_shift_map_figure_(x, y, val, ttl, cb_label, x_label, y_label, cmap, clim_value, FS, marker_size)

fig = figure("Color", "w", "Units", "pixels", "Position", [120 120 720 620]);
ax = axes(fig);

scatter(ax, x, y, marker_size, val, "filled", ...
    "MarkerEdgeAlpha", 0.0, ...
    "MarkerFaceAlpha", 1.0);

axis(ax, "equal");
box(ax, "on");
set(ax, "FontSize", FS, "LineWidth", 1.3, "TickDir", "in", "Box", "on");
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
    delta = max(abs(xmin), 1.0) * 1e-6;
    clim_value = [xmin - delta, xmax + delta];
elseif lower(mode) == "symmetric"
    c = max(abs([xmin, xmax]));
    clim_value = [-c, c];
else
    clim_value = [xmin, xmax];
end

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

function [data, columns, meta] = read_orbital_moment_txt_(txt_file)

fid = fopen(txt_file, "r");
if fid < 0
    error("Cannot open file: %s", txt_file);
end

cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

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
