function plot_chi_vs_doping()

% ============================================================
% USER SETTINGS
% ============================================================
iq_target = 3;
jq_target = -3;

eta = 1e-5;
use_form_factor = true;
boundary_periodic = false;

use_abs_chi = true;   % true: plot |Re chi|, false: plot Re chi

default_dir = 'E:\rg_master\rhombohedral\data';

% ============================================================
% Select folder
% Example:
%   ...\D-0.012\T6.500
% containing:
%   mu0.500000\fermiPatch_xxx.bin
%   mu0.510000\fermiPatch_xxx.bin
% ============================================================
if ~isfolder(default_dir)
    default_dir = pwd;
end

root_dir = uigetdir(default_dir, 'Select T folder containing mu*/fermiPatch*.bin');
if isequal(root_dir, 0)
    return;
end

mu_dirs = dir(fullfile(root_dir, 'mu*'));
mu_dirs = mu_dirs([mu_dirs.isdir]);

if isempty(mu_dirs)
    error('No mu* folders found under: %s', root_dir);
end

doping_list = [];
mu_list = [];
chi_list = [];
chi_imag_list = [];
file_list = strings(0);

fprintf('Selected root:\n%s\n\n', root_dir);
fprintf('q = (%d,%d)\n', iq_target, jq_target);
fprintf('eta = %.3e\n', eta);
fprintf('use_form_factor = %d\n', use_form_factor);
fprintf('boundary_periodic = %d\n\n', boundary_periodic);

% ============================================================
% Loop over mu folders
% ============================================================
for id = 1:numel(mu_dirs)

    mu_dir = fullfile(root_dir, mu_dirs(id).name);

    files = dir(fullfile(mu_dir, 'fermiPatch*.bin'));
    if isempty(files)
        fprintf('[MISS] %s\n', mu_dirs(id).name);
        continue;
    end

    % choose newest bin file in this mu folder
    [~, idx_newest] = max([files.datenum]);
    bin_file = fullfile(mu_dir, files(idx_newest).name);

    fprintf('------------------------------------------------------------\n');
    fprintf('[%d/%d] %s\n', id, numel(mu_dirs), bin_file);

    F = read_fermiPatch_bin_v3(bin_file);

    chi = cal_chi_single_q_from_bin(F, iq_target, jq_target, eta, ...
        use_form_factor, boundary_periodic);

    mu_val = F.EF;
    doping_val = F.doping;

    fprintf('mu     = %.10f\n', mu_val);
    fprintf('doping = %.10f\n', doping_val);
    fprintf('chi_re = %.15e\n', real(chi));
    fprintf('chi_im = %.15e\n', imag(chi));

    mu_list(end+1,1) = mu_val;
    doping_list(end+1,1) = doping_val;
    chi_list(end+1,1) = real(chi);
    chi_imag_list(end+1,1) = imag(chi);
    file_list(end+1,1) = string(bin_file);
end

if isempty(doping_list)
    error('No valid bin files were processed.');
end

% ============================================================
% Sort by doping
% ============================================================
[doping_list, order] = sort(doping_list);
mu_list = mu_list(order);
chi_list = chi_list(order);
chi_imag_list = chi_imag_list(order);
file_list = file_list(order);

if use_abs_chi
    y = abs(chi_list);
    ylab = sprintf('|Re \\chi(q=(%d,%d))|', iq_target, jq_target);
else
    y = chi_list;
    ylab = sprintf('Re \\chi(q=(%d,%d))', iq_target, jq_target);
end

% ============================================================
% Plot
% ============================================================
figure('Color','w');
plot(doping_list, y, '-o', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 6);

xlabel('doping (10^{12} cm^{-2})');
ylabel(ylab);
title(sprintf('\\chi vs doping, q=(%d,%d)', iq_target, jq_target));

set(gca, 'FontSize', 14, 'LineWidth', 1.2);
box on;
grid off;

% ============================================================
% Save MAT result
% ============================================================
OUT.doping = doping_list;
OUT.mu = mu_list;
OUT.chi_real = chi_list;
OUT.chi_imag = chi_imag_list;
OUT.file = file_list;
OUT.iq = iq_target;
OUT.jq = jq_target;
OUT.eta = eta;
OUT.use_form_factor = use_form_factor;
OUT.boundary_periodic = boundary_periodic;

save(fullfile(root_dir, sprintf('chi_vs_doping_q%d_%d.mat', iq_target, jq_target)), 'OUT');

fprintf('\nSaved MAT:\n%s\n', ...
    fullfile(root_dir, sprintf('chi_vs_doping_q%d_%d.mat', iq_target, jq_target)));

end


% ============================================================
% Read fermiPatch bin version 3
% ============================================================
function F = read_fermiPatch_bin_v3(file_path)

fid = fopen(file_path, 'rb');
if fid < 0
    error('Cannot open file: %s', file_path);
end

cleanupObj = onCleanup(@() fclose(fid));

magic   = fread(fid, 1, 'int32');
version = fread(fid, 1, 'int32');
NkTot   = fread(fid, 1, 'int32');
dim     = fread(fid, 1, 'int32');

EF      = fread(fid, 1, 'double');
T_K     = fread(fid, 1, 'double');
doping  = fread(fid, 1, 'double');
filling = fread(fid, 1, 'double');

dx = fread(fid, 1, 'double');
dy = fread(fid, 1, 'double');

mesh_type_raw = fread(fid, 32, '*char')';
mesh_type = erase(char(mesh_type_raw), char(0));
mesh_type = strtrim(mesh_type);

if magic ~= 20260510
    error('Wrong magic number: %d', magic);
end

if version ~= 3
    error('Unsupported version: %d', version);
end

iq = zeros(NkTot, 1);
jq = zeros(NkTot, 1);

kx = zeros(NkTot, 1);
ky = zeros(NkTot, 1);
occ_avg = zeros(NkTot, 1);

evals = zeros(NkTot, dim);
occ_band = zeros(NkTot, dim);

evec_re = zeros(NkTot, dim, dim);
evec_im = zeros(NkTot, dim, dim);

for ik = 1:NkTot
    iq(ik) = fread(fid, 1, 'int32');
    jq(ik) = fread(fid, 1, 'int32');

    kx(ik) = fread(fid, 1, 'double');
    ky(ik) = fread(fid, 1, 'double');
    occ_avg(ik) = fread(fid, 1, 'double');

    evals(ik, :) = fread(fid, dim, 'double')';
    occ_band(ik, :) = fread(fid, dim, 'double')';

    tmp_re = fread(fid, dim * dim, 'double');
    tmp_im = fread(fid, dim * dim, 'double');

    evec_re(ik,:,:) = reshape(tmp_re, [dim, dim]).';
    evec_im(ik,:,:) = reshape(tmp_im, [dim, dim]).';
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
F.evec_re = evec_re;
F.evec_im = evec_im;

F.key_map = containers.Map('KeyType', 'char', 'ValueType', 'double');

for ik = 1:NkTot
    F.key_map(ij_key(iq(ik), jq(ik))) = ik;
end

end


% ============================================================
% Calculate single-q chi
% ============================================================
function chi = cal_chi_single_q_from_bin(F, iq_q, jq_q, eta, ...
    use_form_factor, boundary_periodic)

dim = F.dim;
NkTot = F.NkTot;

chi_sum = complex(0.0, 0.0);
nPair = 0;

[iq_min, iq_max] = bounds(F.iq);
[jq_min, jq_max] = bounds(F.jq);

N1 = iq_max - iq_min + 1;
N2 = jq_max - jq_min + 1;

for ik1 = 1:NkTot

    iq2 = F.iq(ik1) + iq_q;
    jq2 = F.jq(ik1) + jq_q;

    if boundary_periodic && ~strcmp(F.mesh_type, 'hex')
        iq2 = iq_min + mod(iq2 - iq_min, N1);
        jq2 = jq_min + mod(jq2 - jq_min, N2);
    end

    key2 = ij_key(iq2, jq2);

    if ~isKey(F.key_map, key2)
        continue;
    end

    ik2 = F.key_map(key2);

    for b = 1:dim
        Eb = F.evals(ik1, b);
        fb = F.occ_band(ik1, b);

        for m = 1:dim
            Em = F.evals(ik2, m);
            fm = F.occ_band(ik2, m);

            if use_form_factor
                FF = band_form_factor(F, ik1, ik2, b, m);
            else
                FF = 1.0;
            end

            chi_sum = chi_sum + FF * (fm - fb) / complex(Eb - Em, eta);
            nPair = nPair + 1;
        end
    end
end

chi_patch_avg = chi_sum / NkTot;
area_weight = NkTot * F.dx * F.dy;

chi = chi_patch_avg * area_weight;

fprintf('nPair = %d, effective k pairs = %.6f, area_weight = %.6e\n', ...
    nPair, nPair / (dim * dim), area_weight);

end


% ============================================================
% Form factor |<u_kb | u_k+q,m>|^2
% ============================================================
function FF = band_form_factor(F, ik1, ik2, b, m)

u1_re = squeeze(F.evec_re(ik1, :, b));
u1_im = squeeze(F.evec_im(ik1, :, b));

u2_re = squeeze(F.evec_re(ik2, :, m));
u2_im = squeeze(F.evec_im(ik2, :, m));

u1 = complex(u1_re, u1_im);
u2 = complex(u2_re, u2_im);

ov = sum(conj(u1) .* u2);

FF = abs(ov)^2;

end


% ============================================================
% Key for map lookup
% ============================================================
function key = ij_key(iq, jq)

key = sprintf('%d_%d', iq, jq);

end