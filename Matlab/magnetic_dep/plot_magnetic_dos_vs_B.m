function plot_magnetic_dos_vs_B()
clc; close all;

root = fullfile('..', '..', 'data', ...
    'magnetic_dos_sk_60_b1b2_0.00125', ...
    'D-0.010');

B_dirs = list_B_dirs(root);
cmap = hot(256);
B_vals = parse_B_dirs(B_dirs);

figure('Color', 'w', 'Position', [220 220 560 430]);
hold on;

for iB = 1:numel(B_dirs)
    Bdir = B_dirs{iB};
    files = dir(fullfile(root, Bdir, 'magnetic_dos_*.txt'));
    if isempty(files)
        error('No magnetic_dos txt found under %s', fullfile(root, Bdir));
    end

    [~, idx] = max([files.datenum]);
    fpath = fullfile(files(idx).folder, files(idx).name);

    [data, cols, B_T] = read_dos_file(fpath);
    c = interp_color(B_T, B_vals, cmap);

    E = data(:, col_index(cols, 'E_eV'));
    dos_avg = data(:, col_index(cols, 'dos_avg'));
    plot(E, dos_avg, 'LineWidth', 1.2, 'Color', c);
end

box on;
grid on;
set(gca, 'LineWidth', 1.1, 'TickDir', 'in', 'FontSize', 12);
xlabel('energy [eV]');
ylabel('DOS');
title('local DOS vs magnetic field');

colormap(cmap);
clim([min(B_vals), max(B_vals)]);
cb = colorbar;
cb.Label.String = 'B [T]';

savefig(gcf, fullfile(root, 'magnetic_dos_vs_B.fig'));

end

function [data, cols, B_T] = read_dos_file(fpath)
fid = fopen(fpath, 'r');
if fid < 0
    error('Cannot open %s', fpath);
end

cols = {};
B_T = NaN;

while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end

    sline = strtrim(line);
    if startsWith(sline, '# columns =')
        cols = strsplit(strtrim(extractAfter(sline, '=')));
    elseif startsWith(sline, '# B_T =')
        B_T = str2double(strtrim(extractAfter(sline, '=')));
    elseif ~startsWith(sline, '#') && ~isempty(sline)
        fseek(fid, -numel(line)-newline_length(), 'cof');
        break;
    end
end

if isempty(cols)
    error('No column header found in %s', fpath);
end

data = textscan(fid, repmat('%f', 1, numel(cols)), ...
    'CollectOutput', true);
data = data{1};
fclose(fid);
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
