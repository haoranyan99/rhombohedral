function plot_multi_fermi_bin_occ_band()

% ============================================================
% USER SETTINGS
% ============================================================
band_idx = 6;

% "occ", "fs", "energy"
plot_mode = "occ";

marker_size = 35;
FS_tol_eV = 0.01;

a = 2.46;
default_dir = pwd;

% ============================================================
% Manual mu list
% ============================================================

% VHS 1 setting
mu_list = [ ...
    0.75 0.76 ...
    0.772 0.776 0.78 0.784 ...
    0.786 0.788 0.79 0.792 0.794 0.796
];

% % VHS 2 setting
% mu_list = [ ...
%     .56 .57 .58 ...
%     .584 .588 ...
%     .6 .602 .604 .606 .608 ...
%     .61 .62
% ];

mu_list = [0.789 0.79 0.791 0.792];
% mu_list = [0.609 0.61 0.611 0.612];

remove_idx = [];

% Axis limits
xlim_user = [1.35 1.60];
ylim_user = [-0.95 -0.75];

% ============================================================
% Apply removal
% ============================================================
mu_list(remove_idx) = [];
mu_list = sort(mu_list);

supported_n = [1 2 3 4 6 8 9 12];
nFile = numel(mu_list);

if ~ismember(nFile, supported_n)
    error(['Unsupported number of panels = %d\n' ...
           'Use 1,2,3,4,6,8,9,12\n' ...
           'Current mu_list size = %d\n' ...
           'Modify remove_idx manually.'], ...
           nFile, nFile);
end

[nRow, nCol] = layout_for_n_(nFile);

% ============================================================
% Select root directory
% ============================================================
root_dir = uigetdir(default_dir, ...
    'Select root folder containing muXXXX folders');

if isequal(root_dir,0)
    return;
end

% ============================================================
% Load files
% ============================================================
Flist = cell(nFile,1);

for i = 1:nFile

    mu = mu_list(i);
    mu_folder = find_mu_folder_(root_dir, mu);

    if mu_folder == ""
        error('Cannot find folder for mu = %.6f', mu);
    end

    Lbin = dir(fullfile(mu_folder, '*.bin'));

    if isempty(Lbin)
        error('No .bin file found in: %s', mu_folder);
    end

    bin_file = fullfile(Lbin(1).folder, Lbin(1).name);

    fprintf('\n[%d/%d]\n', i, nFile);
    fprintf('target mu = %.6f\n', mu);
    fprintf('folder    = %s\n', mu_folder);
    fprintf('file      = %s\n', bin_file);

    Flist{i} = read_fermiPatch_bin_v3_local_(bin_file);

end

% ============================================================
% Global auto limits
% ============================================================
all_kx = [];
all_ky = [];

for i = 1:nFile
    all_kx = [all_kx; Flist{i}.kx]; %#ok<AGROW>
    all_ky = [all_ky; Flist{i}.ky]; %#ok<AGROW>
end

kmax = max(sqrt(all_kx.^2 + all_ky.^2));

if isempty(xlim_user)
    xlim_use = 1.15 * [-kmax kmax];
else
    xlim_use = xlim_user;
end

if isempty(ylim_user)
    ylim_use = 1.15 * [-kmax kmax];
else
    ylim_use = ylim_user;
end

% ============================================================
% Color settings
% ============================================================
switch lower(plot_mode)

    case 'occ'
        c_range = [0 1];
        cb_label = sprintf('f_{k,%d}', band_idx);

    case 'fs'
        c_range = [0 0.25];
        cb_label = 'f(1-f)';

    case 'energy'
        c_range = [0 FS_tol_eV];
        cb_label = '|E_k - E_F| (eV)';

    otherwise
        error('Unknown plot_mode.');
end

% ============================================================
% Plot
% ============================================================
L = build_RG_lattice_local_(a);

figure( ...
    'Color','w', ...
    'Units','pixels', ...
    'Position',[100 100 420*nCol 380*nRow]);

tiledlayout(nRow, nCol, ...
    'Padding','compact', ...
    'TileSpacing','compact');

for i = 1:nFile

    F = Flist{i};

    occ = F.occ_band(:, band_idx);
    Ek  = F.evals(:, band_idx);

    switch lower(plot_mode)

        case 'occ'
            plot_data = occ;

        case 'fs'
            plot_data = occ .* (1 - occ);

        case 'energy'
            plot_data = abs(Ek - F.EF);
    end

    ax = nexttile;

    hold(ax,'on');
    box(ax,'on');

    % --------------------------------------------------------
    % K-point step coordinate
    %
    % K = (b1 + b2)/3.
    % F.dx and F.dy are fractional mesh step sizes along b1/b2.
    % Therefore K corresponds to approximately:
    % iq_K = round((1/3)/F.dx)
    % jq_K = round((1/3)/F.dy)
    % --------------------------------------------------------
    iq_K = round((1/3) / F.dx);
    jq_K = round((1/3) / F.dy);

    dx_b1_step = F.iq - iq_K;
    dy_b2_step = F.jq - jq_K;

    hSc = scatter(ax, ...
        F.kx, ...
        F.ky, ...
        marker_size, ...
        plot_data, ...
        'filled', ...
        'MarkerEdgeColor','none');

    % --------------------------------------------------------
    % Minimal custom datatip
    % --------------------------------------------------------
    
    rows = hSc.DataTipTemplate.DataTipRows;
    
    % overwrite rows
    rows(1) = dataTipTextRow('dx_b1_step', dx_b1_step);
    rows(2) = dataTipTextRow('dy_b2_step', dy_b2_step);
    
    % append extra rows
    rows(end+1) = dataTipTextRow('energy', Ek);
    rows(end+1) = dataTipTextRow('occ', occ);
    
    hSc.DataTipTemplate.DataTipRows = rows;

    plot_BZ_boundary_local_(ax, L);

    % --------------------------------------------------------
    % Draw b1 and b2 directions
    % --------------------------------------------------------
    
    origin = [1.35, -0.75];   % choose a nice corner
    
    scale_arrow = 0.1;
    
    b1v = scale_arrow * L.b1 / norm(L.b1);
    b2v = scale_arrow * L.b2 / norm(L.b2);
    
    quiver(ax, ...
        origin(1), origin(2), ...
        b1v(1), b1v(2), ...
        0, ...
        'w', ...
        'LineWidth', 2, ...
        'MaxHeadSize', 1.5);
    
    quiver(ax, ...
        origin(1), origin(2), ...
        b2v(1), b2v(2), ...
        0, ...
        'c', ...
        'LineWidth', 2, ...
        'MaxHeadSize', 1.5);
    
    text(ax, ...
        origin(1)+b1v(1)*1.1, ...
        origin(2)+b1v(2)*1.1, ...
        'b_1', ...
        'Color','w', ...
        'FontSize',12, ...
        'FontWeight','bold');
    
    text(ax, ...
        origin(1)+b2v(1)*1.1, ...
        origin(2)+b2v(2)*1.1, ...
        'b_2', ...
        'Color','c', ...
        'FontSize',12, ...
        'FontWeight','bold');

    axis(ax,'equal');

    xlim(ax, xlim_use);
    ylim(ax, ylim_use);

    caxis(ax, c_range);
    colormap(ax, flip(hot));

    title(ax, ...
        sprintf('\\mu=%.4f, dop=%.4f', ...
        F.EF, F.doping), ...
        'FontSize', 12, ...
        'FontWeight','normal');

    xlabel(ax, 'k_x');
    ylabel(ax, 'k_y');

    set(ax, ...
        'FontSize', 11, ...
        'LineWidth', 1.2, ...
        'TickDir','in');

end

cb = colorbar;
cb.Layout.Tile = 'east';

ylabel(cb, cb_label, ...
    'FontSize', 15);

sgtitle( ...
    sprintf('Band %d, mode = %s', ...
    band_idx, plot_mode), ...
    'FontSize', 18, ...
    'FontWeight','bold');

end


% ============================================================
% Layout helper
% ============================================================
function [nRow, nCol] = layout_for_n_(n)

switch n

    case 1
        nRow = 1; nCol = 1;

    case 2
        nRow = 1; nCol = 2;

    case 3
        nRow = 1; nCol = 3;

    case 4
        nRow = 2; nCol = 2;

    case 6
        nRow = 2; nCol = 3;

    case 8
        nRow = 2; nCol = 4;

    case 9
        nRow = 3; nCol = 3;

    case 12
        nRow = 3; nCol = 4;

    otherwise
        error('Unsupported number.');

end

end


% ============================================================
% Find mu folder
% ============================================================
function mu_folder = find_mu_folder_(root_dir, mu)

mu_folder = "";

L = dir(fullfile(root_dir, 'mu*'));

best_err = inf;

for i = 1:numel(L)

    if ~L(i).isdir
        continue;
    end

    name = string(L(i).name);

    tok = regexp(name, ...
        '^mu([+\-]?\d*\.?\d+)$', ...
        'tokens', 'once');

    if isempty(tok)
        continue;
    end

    mu_here = str2double(tok{1});

    err = abs(mu_here - mu);

    if err < best_err

        best_err = err;

        mu_folder = string(fullfile( ...
            L(i).folder, ...
            L(i).name));
    end

end

if best_err > 1e-6
    mu_folder = "";
end

end


% ============================================================
% Read fermiPatch bin
% ============================================================
function F = read_fermiPatch_bin_v3_local_(file_path)

fid = fopen(file_path, 'rb');

if fid < 0
    error('Cannot open file: %s', file_path);
end

cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

magic   = fread(fid, 1, 'int32');
version = fread(fid, 1, 'int32');

NkTot = fread(fid, 1, 'int32');
dim   = fread(fid, 1, 'int32');

EF      = fread(fid, 1, 'double');
T_K     = fread(fid, 1, 'double');
doping  = fread(fid, 1, 'double');
filling = fread(fid, 1, 'double');

dx = fread(fid, 1, 'double');
dy = fread(fid, 1, 'double');

mesh_raw = fread(fid, 32, '*char')';
mesh_type = char(mesh_raw);
mesh_type = erase(mesh_type, char(0));
mesh_type = strtrim(mesh_type);

iq = zeros(NkTot,1);
jq = zeros(NkTot,1);
kx = zeros(NkTot,1);
ky = zeros(NkTot,1);
occ_avg = zeros(NkTot,1);
evals = zeros(NkTot, dim);
occ_band = zeros(NkTot, dim);

for ik = 1:NkTot

    iq(ik) = fread(fid,1,'int32');
    jq(ik) = fread(fid,1,'int32');

    kx(ik) = fread(fid,1,'double');
    ky(ik) = fread(fid,1,'double');

    occ_avg(ik) = fread(fid,1,'double');

    evals(ik,:) = fread(fid, dim, 'double')';
    occ_band(ik,:) = fread(fid, dim, 'double')';

    fread(fid, dim*dim, 'double');
    fread(fid, dim*dim, 'double');

end

F.magic = magic;
F.version = version;
F.NkTot = NkTot;
F.dim = dim;
F.EF = EF;
F.T_K = T_K;
F.doping = doping;
F.filling = filling;
F.dx = dx;
F.dy = dy;
F.mesh_type = mesh_type;
F.iq = iq;
F.jq = jq;
F.kx = kx;
F.ky = ky;
F.occ_avg = occ_avg;
F.evals = evals;
F.occ_band = occ_band;

end


% ============================================================
% Build reciprocal lattice
% ============================================================
function L = build_RG_lattice_local_(a)

a1 = [0.0, -a];
a2 = [sqrt(3)/2 * a, 0.5 * a];

area = a1(1)*a2(2) - a1(2)*a2(1);

b1 = 2*pi * [ a2(2), -a2(1)] / area;
b2 = 2*pi * [-a1(2),  a1(1)] / area;

L.b1 = b1;
L.b2 = b2;

end


% ============================================================
% Plot first BZ boundary
% ============================================================
function plot_BZ_boundary_local_(ax, L)

b1 = L.b1;
b2 = L.b2;

corners = [
     ( b1 + b2)/3;
     (-b1 + 2*b2)/3;
     (-2*b1 + b2)/3;
    -( b1 + b2)/3;
      (b1 - 2*b2)/3;
      (2*b1 - b2)/3;
      (b1 + b2)/3
];

plot(ax, corners(:,1), corners(:,2), ...
    'c-', 'LineWidth', 1.8);

end