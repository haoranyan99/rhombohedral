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

% Scan mode:
%   "mu"     : plot different mu values for a selected branch/group.
%   "branch" : plot different branches at one selected mu.
scan_mode = "branch";
fixed_mu = 0.791;

% Branch labels:
%   1 = valley_plus  spin_up
%   2 = valley_plus  spin_down
%   3 = valley_minus spin_up
%   4 = valley_minus spin_down
%
% In scan_mode="mu", multiple entries are averaged for every mu.
% In scan_mode="branch", each entry becomes one panel.
branch_list = [1 2 3 4];

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
% [] means auto-fit to the selected fermi patch.
xlim_user = [];
ylim_user = [];

% ============================================================
% Apply removal
% ============================================================
mu_list(remove_idx) = [];
mu_list = sort(mu_list);

[panel_mu_list, panel_branch_groups, panel_labels] = build_panel_cases_( ...
    scan_mode, mu_list, fixed_mu, branch_list);

supported_n = [1 2 3 4 6 8 9 12];
nFile = numel(panel_mu_list);

if ~ismember(nFile, supported_n)
    error(['Unsupported number of panels = %d\n' ...
           'Use 1,2,3,4,6,8,9,12\n' ...
           'Current panel count = %d\n' ...
           'Modify mu_list/remove_idx or branch_list.'], ...
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

    mu = panel_mu_list(i);
    branch_group = panel_branch_groups{i};
    mu_folder = find_mu_folder_(root_dir, mu);

    if mu_folder == ""
        error('Cannot find folder for mu = %.6f', mu);
    end

    bin_files = find_branch_bin_files_(mu_folder, branch_group);

    fprintf('\n[%d/%d]\n', i, nFile);
    fprintf('target mu = %.6f\n', mu);
    fprintf('panel     = %s\n', panel_labels(i));
    fprintf('folder    = %s\n', mu_folder);
    fprintf('branches  = %d\n', numel(bin_files));
    for ib = 1:numel(bin_files)
        fprintf('file      = %s\n', bin_files(ib));
    end

    Flist{i} = read_and_average_fermi_files_(bin_files);
    Flist{i}.panel_label = panel_labels(i);

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

all_kx = all_kx(isfinite(all_kx));
all_ky = all_ky(isfinite(all_ky));

if isempty(all_kx) || isempty(all_ky)
    error('No finite k points found. Check bin reader and selected files.');
end

pad_x = 0.06 * max(eps, max(all_kx) - min(all_kx));
pad_y = 0.06 * max(eps, max(all_ky) - min(all_ky));

if isempty(xlim_user)
    xlim_use = [min(all_kx) - pad_x, max(all_kx) + pad_x];
else
    xlim_use = xlim_user;
end

if isempty(ylim_user)
    ylim_use = [min(all_ky) - pad_y, max(all_ky) + pad_y];
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

    if band_idx > F.dim
        error('Requested band_idx=%d but bin file has dim=%d.', band_idx, F.dim);
    end

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

    ok_plot = isfinite(F.kx) & isfinite(F.ky) & isfinite(plot_data) & ...
              isfinite(Ek) & isfinite(occ);
    if ~any(ok_plot)
        error('No finite k/plot data for panel %d. Check bin reader/version.', i);
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
        F.kx(ok_plot), ...
        F.ky(ok_plot), ...
        marker_size, ...
        plot_data(ok_plot), ...
        'filled', ...
        'MarkerEdgeColor','none');

    % --------------------------------------------------------
    % Minimal custom datatip
    % --------------------------------------------------------
    
    rows = hSc.DataTipTemplate.DataTipRows;
    
    % overwrite rows
    rows(1) = dataTipTextRow('dx_b1_step', dx_b1_step(ok_plot));
    rows(2) = dataTipTextRow('dy_b2_step', dy_b2_step(ok_plot));
    
    % append extra rows
    rows(end+1) = dataTipTextRow('energy', Ek(ok_plot));
    rows(end+1) = dataTipTextRow('occ', occ(ok_plot));
    
    hSc.DataTipTemplate.DataTipRows = rows;

    plot_BZ_boundary_local_(ax, L);

    % --------------------------------------------------------
    % Draw b1 and b2 directions
    % --------------------------------------------------------
    
    origin = [ ...
        xlim_use(1) + 0.07 * diff(xlim_use), ...
        ylim_use(1) + 0.10 * diff(ylim_use)];

    scale_arrow = 0.12 * min(diff(xlim_use), diff(ylim_use));
    
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

    xlim(ax, xlim_use);
    ylim(ax, ylim_use);
    axis(ax,'equal');

    caxis(ax, c_range);
    colormap(ax, flip(hot));

    title(ax, ...
        sprintf('%s, dop=%.4f', ...
        F.panel_label, F.doping), ...
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
    sprintf('Band %d, mode = %s, scan = %s', ...
    band_idx, plot_mode, scan_mode), ...
    'FontSize', 18, ...
    'FontWeight','bold');

end


% ============================================================
% Layout helper
% ============================================================
function [panel_mu_list, panel_branch_groups, panel_labels] = build_panel_cases_( ...
    scan_mode, mu_list, fixed_mu, branch_list)

scan_mode = lower(string(scan_mode));
branch_list = branch_list(:).';
validate_branch_list_(branch_list);

switch scan_mode
    case "mu"
        panel_mu_list = mu_list(:).';
        panel_branch_groups = cell(numel(panel_mu_list), 1);
        panel_labels = strings(numel(panel_mu_list), 1);
        branch_label = branch_group_label_(branch_list);

        for i = 1:numel(panel_mu_list)
            panel_branch_groups{i} = branch_list;
            panel_labels(i) = sprintf("\\mu=%.4f, %s", ...
                panel_mu_list(i), branch_label);
        end

    case "branch"
        panel_mu_list = fixed_mu * ones(1, numel(branch_list));
        panel_branch_groups = cell(numel(branch_list), 1);
        panel_labels = strings(numel(branch_list), 1);

        for i = 1:numel(branch_list)
            panel_branch_groups{i} = branch_list(i);
            panel_labels(i) = sprintf("\\mu=%.4f, %s", ...
                fixed_mu, branch_label_(branch_list(i)));
        end

    otherwise
        error('Unknown scan_mode=%s. Use "mu" or "branch".', scan_mode);
end

end

function validate_branch_list_(branch_list)

if isempty(branch_list) || any(~ismember(branch_list, 1:4))
    error('branch_list must contain branch numbers from 1 to 4.');
end

end

function label = branch_group_label_(branch_list)

if numel(branch_list) == 1
    label = branch_label_(branch_list(1));
else
    parts = strings(1, numel(branch_list));
    for i = 1:numel(branch_list)
        parts(i) = "b" + string(branch_list(i));
    end
    label = "avg(" + strjoin(parts, "+") + ")";
end

end

function label = branch_label_(branch_id)

switch branch_id
    case 1
        label = "b1: K+, up";
    case 2
        label = "b2: K+, down";
    case 3
        label = "b3: K-, up";
    case 4
        label = "b4: K-, down";
    otherwise
        error('Unknown branch_id=%d.', branch_id);
end

end

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
% Branch file selection and averaging
% ============================================================
function bin_files = find_branch_bin_files_(mu_folder, branch_group)

bin_files = strings(0,1);
branch_group = branch_group(:).';

for ib = 1:numel(branch_group)
    pattern = branch_file_pattern_(branch_group(ib));
    L = dir(fullfile(mu_folder, pattern));
    if isempty(L)
        error('Missing branch file in %s: %s', mu_folder, pattern);
    end
    bin_files(end+1,1) = string(fullfile(L(1).folder, L(1).name)); %#ok<AGROW>
end

end

function pattern = branch_file_pattern_(branch_id)

switch branch_id
    case 1
        pattern = 'fermi_valley_plus_spin_up_patch.bin';
    case 2
        pattern = 'fermi_valley_plus_spin_down_patch.bin';
    case 3
        pattern = 'fermi_valley_minus_spin_up_patch.bin';
    case 4
        pattern = 'fermi_valley_minus_spin_down_patch.bin';
    otherwise
        error('Unknown branch_id=%d.', branch_id);
end

end

function Favg = read_and_average_fermi_files_(files)

files = string(files(:));
Favg = read_fermiPatch_bin_v3_local_(files(1));
Favg.source_files = files;
Favg.n_averaged = numel(files);

if numel(files) == 1
    return;
end

eval_sum = Favg.evals;
occ_band_sum = Favg.occ_band;
occ_avg_sum = Favg.occ_avg;
doping_sum = Favg.doping;
filling_sum = Favg.filling;
doping_spin_sum = Favg.doping_spin;
filling_spin_sum = Favg.filling_spin;

for i = 2:numel(files)
    F = read_fermiPatch_bin_v3_local_(files(i));
    assert_same_kmesh_(Favg, F, files(i));

    eval_sum = eval_sum + F.evals;
    occ_band_sum = occ_band_sum + F.occ_band;
    occ_avg_sum = occ_avg_sum + F.occ_avg;
    doping_sum = doping_sum + F.doping;
    filling_sum = filling_sum + F.filling;
    doping_spin_sum = doping_spin_sum + F.doping_spin;
    filling_spin_sum = filling_spin_sum + F.filling_spin;
end

n = numel(files);
Favg.evals = eval_sum ./ n;
Favg.occ_band = occ_band_sum ./ n;
Favg.occ_avg = occ_avg_sum ./ n;
Favg.doping = doping_sum ./ n;
Favg.filling = filling_sum ./ n;
Favg.doping_spin = doping_spin_sum ./ n;
Favg.filling_spin = filling_spin_sum ./ n;
Favg.spin_sign = 0;

end

function assert_same_kmesh_(A, B, file_path)

if A.NkTot ~= B.NkTot || A.dim ~= B.dim
    error('Cannot average branch with different size/dim: %s', file_path);
end

tol = 1e-10;
same_grid = isequal(A.iq, B.iq) && isequal(A.jq, B.jq) && ...
            max(abs(A.kx - B.kx)) < tol && ...
            max(abs(A.ky - B.ky)) < tol;

if ~same_grid
    error('Cannot average branch with different k mesh: %s', file_path);
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
if magic ~= 20260510
    error('Unexpected fermiPatch magic=%d in %s.', magic, file_path);
end

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

spin_sign = 0;
doping_spin = NaN;
filling_spin = NaN;
if version >= 5
    spin_sign = fread(fid, 1, 'int32');
    doping_spin = fread(fid, 1, 'double');
    filling_spin = fread(fid, 1, 'double');
end

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
F.spin_sign = spin_sign;
F.doping_spin = doping_spin;
F.filling_spin = filling_spin;
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
