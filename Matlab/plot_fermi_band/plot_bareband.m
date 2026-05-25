function plot_bareband()

clc; close all;

% =============================
% Choose file
% =============================
default_root = "/Users/haoranyan/rg_master/data/";
if ~isfolder(default_root)
    default_root = pwd;
end

[fname, fpath] = uigetfile({'*.txt','Band data (*.txt)'}, ...
    'Select RG band file', fullfile(default_root, '*.txt'));

if isequal(fname,0)
    return;
end

in_path = fullfile(fpath, fname);

% =============================
% Read data
% =============================
A = readmatrix(in_path, 'CommentStyle', '#');

if isempty(A) || size(A,2) < 5
    error('Invalid data.');
end

s = A(:,2);
E = A(:,5:end);
norb = size(E,2);

fprintf("Detected %d bands.\n", norb);

% =============================
% Shift mu: middle of two central flat bands E5 and E6
% File columns are E0, E1, ..., so:
% E5 -> E(:,6)
% E6 -> E(:,7)
% =============================
if norb < 7
    error('Need at least E6, i.e. at least 7 bands: E0...E6.');
end

E5 = E(:,6);
E6 = E(:,7);

E5_max = max(E5, [], 'omitnan');
E6_min = min(E6, [], 'omitnan');

mu_shift = 0.5 * (E5_max + E6_min);
E = E - mu_shift;

fprintf("mu shift = %.10f eV\n", mu_shift);
fprintf("After shift: max(E5)=%.10f eV, min(E6)=%.10f eV\n", ...
    E5_max - mu_shift, E6_min - mu_shift);

% =============================
% Parse xticks
% =============================
xtick_pos = [];
xtick_lab = {};

fid = fopen(in_path, 'r');
while ~feof(fid)
    line = strtrim(fgetl(fid));

    if startsWith(line, '# xticks:')
        blocks = regexp(line, '\(([^)]+)\)', 'tokens');

        for i = 1:numel(blocks)
            parts = strsplit(strtrim(blocks{i}{1}), ',');

            if numel(parts) >= 2
                xtick_lab{end+1} = strtrim(parts{1});
                xtick_pos(end+1) = str2double(parts{2});
            end
        end
        break;
    end
end
fclose(fid);

if isempty(xtick_pos)
    xtick_pos = linspace(s(1), s(end), 3);
    xtick_lab = {'start','mid','end'};
end

% 强制边界对齐
xtick_pos(1)   = s(1);
xtick_pos(end) = s(end);

% =============================
% Plot
% =============================
figure('Color','w','Position',[120 120 900 520]);
hold on; box on; grid on;

for ib = 1:norb
    plot(s, E(:,ib), 'LineWidth', 1.2);
end

xlabel('k-line coordinate s (1/\AA)');
ylabel('Energy (eV)');
set(gca, ...
    'FontSize', 13, ...
    'LineWidth', 1.1, ...
    'XTick', xtick_pos, ...
    'XTickLabel', xtick_lab, ...
    'XLim', [s(1), s(end)]);

% 高对称点竖线
for i = 2:numel(xtick_pos)-1
    xline(xtick_pos(i), '-', 'Color', [0.75 0.75 0.75], ...
        'LineWidth', 0.8, 'HandleVisibility','off');
end

end