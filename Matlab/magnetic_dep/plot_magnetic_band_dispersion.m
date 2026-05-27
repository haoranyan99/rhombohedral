function plot_magnetic_band_dispersion()

root = fullfile('..', '..', 'data', ...
    'magnetic_band_sk_localK_MKKp_Nkseg1000_frac0.1_D-0.01', ...
    'D-0.010');

B_dirs = list_B_dirs(root);
band_list = [4 5 6];
energy_window_eV = [0.45 0.85];

spin_data = collect_spin_band_data(root, B_dirs, band_list);
cmap = hot(256);

plot_one_spin(root, spin_data, 'up', cmap, energy_window_eV);
plot_one_spin(root, spin_data, 'down', cmap, energy_window_eV);

end

function S = collect_spin_band_data(root, B_dirs, band_list)
S = struct();
S.B_vals = parse_B_dirs(B_dirs);
S.entries = struct([]);
S.xtick_pos = [];
S.xtick_lab = {};

for iB = 1:numel(B_dirs)
    Bdir = B_dirs{iB};
    files = dir(fullfile(root, Bdir, 'magnetic_band_*.txt'));
    if isempty(files)
        error('No magnetic_band txt found under %s', fullfile(root, Bdir));
    end

    [~, idx] = max([files.datenum]);
    fpath = fullfile(files(idx).folder, files(idx).name);

    [data, cols, xtick_pos, xtick_lab, B_T] = read_band_file(fpath);

    entry = struct();
    entry.B_T = B_T;
    entry.s = data(:, col_index(cols, 's'));
    entry.Eup = cell(size(band_list));
    entry.Edn = cell(size(band_list));

    for ib = 1:numel(band_list)
        b = band_list(ib);
        entry.Eup{ib} = data(:, col_index(cols, sprintf('Eup_b%d', b)));
        entry.Edn{ib} = data(:, col_index(cols, sprintf('Edn_b%d', b)));
    end

    S.entries(end + 1) = entry; %#ok<AGROW>
    S.xtick_pos = xtick_pos;
    S.xtick_lab = xtick_lab;
end
end

function plot_one_spin(root, S, spin_name, cmap, energy_window_eV)
fig = figure('Color', 'w', 'Position', [200 200 720 520]);
ax = axes(fig);
hold(ax, 'on');

for iB = 1:numel(S.entries)
    entry = S.entries(iB);
    c = interp_color(entry.B_T, S.B_vals, cmap);

    if strcmp(spin_name, 'up')
        bands = entry.Eup;
    else
        bands = entry.Edn;
    end

    for ib = 1:numel(bands)
        plot(ax, entry.s, bands{ib}, '-', 'LineWidth', 0.9, ...
            'Color', c, 'HandleVisibility', 'off');
    end
end

first_s = S.entries(1).s;
xlim(ax, [min(first_s), max(first_s)]);
ylim(ax, energy_window_eV);
grid(ax, 'on');
box(ax, 'on');
set(ax, 'LineWidth', 1.1, 'TickDir', 'in', 'FontSize', 12);

if ~isempty(S.xtick_pos)
    xticks(ax, S.xtick_pos);
    xticklabels(ax, S.xtick_lab);
end

colormap(fig, cmap);
caxis(ax, [min(S.B_vals), max(S.B_vals)]);
cb = colorbar(ax);
cb.Label.String = 'B [T]';

xlabel(ax, 'k path');
ylabel(ax, 'energy [eV]');
title(ax, sprintf('spin %s', spin_name));

savefig(fig, fullfile(root, sprintf('magnetic_band_dispersion_spin_%s.fig', spin_name)));
end

function c = interp_color(B_T, B_vals, cmap)
B_min = min(B_vals);
B_max = max(B_vals);
if abs(B_max - B_min) < eps
    idx = 1;
else
    idx = 1 + round((size(cmap, 1) - 1) * (B_T - B_min) / (B_max - B_min));
end
idx = max(1, min(size(cmap, 1), idx));
c = cmap(idx, :);
end

function B_dirs = list_B_dirs(root)
d = dir(fullfile(root, 'B*T'));
d = d([d.isdir]);
B = zeros(size(d));
for i = 1:numel(d)
    B(i) = parse_B(d(i).name);
end
[~, order] = sort(B);
B_dirs = {d(order).name};
end

function B = parse_B_dirs(B_dirs)
B = zeros(1, numel(B_dirs));
for i = 1:numel(B_dirs)
    B(i) = parse_B(B_dirs{i});
end
end

function B = parse_B(name)
s = erase(erase(name, 'B'), 'T');
s = strrep(s, 'p', '.');
B = str2double(s);
end

function [data, cols, xtick_pos, xtick_lab, B_T] = read_band_file(fpath)
fid = fopen(fpath, 'r');
if fid < 0
    error('Cannot open %s', fpath);
end

cols = {};
xtick_pos = [];
xtick_lab = {};
B_T = NaN;

while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end

    if startsWith(line, '# columns =')
        cols = strsplit(strtrim(extractAfter(line, '=')));
    elseif startsWith(line, '# xtick =')
        tok = strsplit(strtrim(extractAfter(line, '=')));
        xtick_pos(end + 1) = str2double(tok{1}); %#ok<AGROW>
        xtick_lab{end + 1} = tok{2}; %#ok<AGROW>
    elseif startsWith(line, '# B_T =')
        B_T = str2double(strtrim(extractAfter(line, '=')));
    elseif ~startsWith(strtrim(line), '#') && ~isempty(strtrim(line))
        fseek(fid, -numel(line)-newline_length(), 'cof');
        break;
    end
end

data = textscan(fid, repmat('%f', 1, numel(cols)), ...
    'CollectOutput', true);
data = data{1};
fclose(fid);
end

function idx = col_index(cols, name)
idx = find(strcmp(cols, name), 1);
if isempty(idx)
    error('Column not found: %s', name);
end
end

function n = newline_length()
if ispc
    n = 2;
else
    n = 1;
end
end
