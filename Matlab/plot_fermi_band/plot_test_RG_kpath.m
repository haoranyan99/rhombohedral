% File: plot_test_RG_kpath.m
% Plot k-space path Γ–M–K–Γ on top of primitive reciprocal cell (parallelogram)

clear; clc;

fn_bz   = 'data/test_RG_kpath_bz_cell.txt';
fn_path = 'data/test_RG_kpath_GMKG.txt';
fn_hs   = 'data/test_RG_kpath_highsym.txt';

% -------- load BZ cell vertices (vx vy) --------
bz = readmatrix(fn_bz, 'CommentStyle', '#');  % expected: N x 2 (closed polygon, last point = first)

% -------- load path (ik kx ky s) --------
path = readmatrix(fn_path, 'CommentStyle', '#');
kx = path(:,2);
ky = path(:,3);

% -------- load high symmetry points (name kx ky s) --------
fid = fopen(fn_hs,'r');
if fid < 0
    error('Cannot open %s', fn_hs);
end

names = {};
kx_h  = [];
ky_h  = [];
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line), continue; end
    if startsWith(line, '#'), continue; end

    tok = strsplit(line);
    if numel(tok) < 3
        continue;
    end

    names{end+1} = tok{1}; %#ok<SAGROW>
    kx_h(end+1)  = str2double(tok{2}); %#ok<SAGROW>
    ky_h(end+1)  = str2double(tok{3}); %#ok<SAGROW>
end
fclose(fid);

% Map "G" -> "Γ" for nicer labels
for i = 1:numel(names)
    if strcmp(names{i}, 'G')
        names{i} = char(915); % Γ
    end
end

% -------- plot --------
figure('Color','w'); hold on; axis equal;

% primitive reciprocal cell
plot(bz(:,1), bz(:,2), 'k--', 'LineWidth', 1.2);

% k-path
plot(kx, ky, 'r-', 'LineWidth', 2.0);

% high-symmetry markers
scatter(kx_h, ky_h, 60, 'k', 'filled');

% labels with small offset based on axis span
xmin = min([bz(:,1); kx]); xmax = max([bz(:,1); kx]);
ymin = min([bz(:,2); ky]); ymax = max([bz(:,2); ky]);
dx = 0.02 * (xmax - xmin);
dy = 0.02 * (ymax - ymin);

for i = 1:numel(names)
    text(kx_h(i) + dx, ky_h(i) + dy, names{i}, ...
        'FontSize', 14, 'FontWeight', 'bold');
end

xlabel('k_x (1/Å)');
ylabel('k_y (1/Å)');
title('\Gamma–M–K–\Gamma path on reciprocal primitive cell');
grid on;

% nice limits with padding
padx = 0.08 * (xmax - xmin);
pady = 0.08 * (ymax - ymin);
xlim([xmin - padx, xmax + padx]);
ylim([ymin - pady, ymax + pady]);

mkdir('plot');
saveas(gcf, 'plot/test_RG_kpath.png');
