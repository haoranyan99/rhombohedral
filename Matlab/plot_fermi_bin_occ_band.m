function plot_fermi_bin_occ_band()

clc;
close all;

% ============================================================
% USER SETTINGS
% ============================================================
band_idx = 6;

% "occ"    -> occupation
% "fs"     -> f(1-f)
% "energy" -> |E-EF|
plot_mode = "fs";

a = 2.46;

default_dir = pwd;

marker_size = 45;

FS_tol_eV = 0.01;

% ============================================================
% Select bin
% ============================================================
[fname, fpath] = uigetfile( ...
    fullfile(default_dir, '*.bin'), ...
    'Select fermiPatch_*.bin');

if isequal(fname,0)
    return;
end

bin_file = fullfile(fpath, fname);

% ============================================================
% Read bin
% ============================================================
F = read_fermiPatch_bin_v3(bin_file);

fprintf('Read: %s\n', bin_file);
fprintf('NkTot     = %d\n', F.NkTot);
fprintf('dim       = %d\n', F.dim);
fprintf('EF        = %.10f eV\n', F.EF);
fprintf('T_K       = %.6f K\n', F.T_K);
fprintf('doping    = %.10f\n', F.doping);
fprintf('filling   = %.10f\n', F.filling);
fprintf('mesh_type = %s\n', F.mesh_type);

if band_idx < 1 || band_idx > F.dim
    error('band_idx must be between 1 and %d', F.dim);
end

% ============================================================
% Select plotted quantity
% ============================================================
occ = F.occ_band(:, band_idx);
Ek  = F.evals(:, band_idx);

switch lower(plot_mode)

    case 'occ'

        plot_data = occ;
        cb_label = sprintf('f_{k,%d}', band_idx);
        c_range = [0 1];

    case 'fs'

        plot_data = occ .* (1 - occ);

        cb_label = 'f(1-f)';

        c_range = [0 0.25];

    case 'energy'

        plot_data = abs(Ek - F.EF);

        cb_label = '|E_k - E_F| (eV)';

        c_range = [0 FS_tol_eV];

    otherwise

        error('Unknown plot_mode');

end

% ============================================================
% Build lattice
% ============================================================
L = build_RG_lattice(a);

% ============================================================
% Plot
% ============================================================
figure('Color','w');

hold on;
box on;

scatter( ...
    F.kx, ...
    F.ky, ...
    marker_size, ...
    plot_data, ...
    'filled', ...
    'MarkerEdgeColor', 'none');

% ============================================================
% Colormap
% ============================================================
colormap(flip(hot));

caxis(c_range);

cb = colorbar;

ylabel(cb, cb_label, ...
    'FontSize', 16);

% ============================================================
% BZ
% ============================================================
plot_BZ_boundary(L);

% ============================================================
% High symmetry points
% ============================================================
scatter(L.Gamma(1), L.Gamma(2), ...
    80, 'c', 'filled');

text( ...
    L.Gamma(1), ...
    L.Gamma(2), ...
    '  \Gamma', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

scatter(L.M(1), L.M(2), ...
    80, 'c', 'filled');

text( ...
    L.M(1), ...
    L.M(2), ...
    '  M', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

scatter(L.K(1), L.K(2), ...
    80, 'c', 'filled');

text( ...
    L.K(1), ...
    L.K(2), ...
    '  K', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

scatter(L.Kp(1), L.Kp(2), ...
    80, 'c', 'filled');

text( ...
    L.Kp(1), ...
    L.Kp(2), ...
    '  K''', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

% ============================================================
% Axis
% ============================================================
axis equal;

xlabel('k_x', 'FontSize', 18);
ylabel('k_y', 'FontSize', 18);

title( ...
    sprintf( ...
        'Band %d   EF = %.4f eV', ...
        band_idx, ...
        F.EF), ...
    'FontSize', 18);

set(gca, ...
    'FontSize', 15, ...
    'LineWidth', 1.8, ...
    'TickDir', 'in');

% ============================================================
% Auto limits
% ============================================================
kmax = max( ...
    sqrt(F.kx.^2 + F.ky.^2));

xlim(1.15 * [-kmax, kmax]);
ylim(1.15 * [-kmax, kmax]);

end


% ============================================================
% Read bin
% ============================================================
function F = read_fermiPatch_bin_v3(file_path)

fid = fopen(file_path, 'rb');

if fid < 0
    error('Cannot open file');
end

cleanupObj = onCleanup(@() fclose(fid));

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
% RG lattice
% ============================================================
function L = build_RG_lattice(a)

a1 = [0.0, -a];
a2 = [sqrt(3)/2 * a, 0.5 * a];

area = ...
    a1(1)*a2(2) ...
  - a1(2)*a2(1);

twoPi = 2*pi;

b1 = twoPi * [ a2(2), -a2(1)] / area;
b2 = twoPi * [-a1(2),  a1(1)] / area;

Gamma = [0 0];

M  = 0.5 * b1;
K  = (b1 + b2)/3;
Kp = 2*(b1 + b2)/3;

L.a1 = a1;
L.a2 = a2;

L.b1 = b1;
L.b2 = b2;

L.Gamma = Gamma;
L.M  = M;
L.K  = K;
L.Kp = Kp;

end


% ============================================================
% Plot BZ
% ============================================================
function plot_BZ_boundary(L)

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

plot( ...
    corners(:,1), ...
    corners(:,2), ...
    'c-', ...
    'LineWidth', 2.5);

% ============================================================
% b1
% ============================================================
quiver( ...
    0,0, ...
    b1(1), ...
    b1(2), ...
    0, ...
    'Color', [1 0 0], ...
    'LineWidth', 2.5, ...
    'MaxHeadSize', 0.5);

text( ...
    1.05*b1(1), ...
    1.05*b1(2), ...
    'b_1', ...
    'Color', [1 0 0], ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

% ============================================================
% b2
% ============================================================
quiver( ...
    0,0, ...
    b2(1), ...
    b2(2), ...
    0, ...
    'Color', [0 0.3 1], ...
    'LineWidth', 2.5, ...
    'MaxHeadSize', 0.5);

text( ...
    1.05*b2(1), ...
    1.05*b2(2), ...
    'b_2', ...
    'Color', [0 0.3 1], ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

end