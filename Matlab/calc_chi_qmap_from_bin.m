function calc_chi_qmap_from_bin()

% ============================================================
% USER SETTINGS
% qx,qy are in units of b1,b2:
%   q = qx * b1 + qy * b2
%
% bin grid:
%   qx = iq * F.dx
%   qy = jq * F.dy
% ============================================================

qx_min = -0.005;
qx_max =  0;
qy_min =  0;
qy_max =  0.005;

eta = 1e-5;
use_form_factor = true;
boundary_periodic = false;

plot_quantity = "abs";   % "real", "abs", "imag"

% ============================================================
% Select bin file
% ============================================================

[fname, fpath] = uigetfile('*.bin', 'Select fermiPatch bin file');
if isequal(fname, 0)
    return;
end

bin_file = fullfile(fpath, fname);

% ============================================================
% Read bin
% ============================================================

F = read_fermiPatch_bin_v3(bin_file);

fprintf('Read file: %s\n', bin_file);
fprintf('NkTot     = %d\n', F.NkTot);
fprintf('dim       = %d\n', F.dim);
fprintf('EF        = %.10f eV\n', F.EF);
fprintf('T_K       = %.6f K\n', F.T_K);
fprintf('doping    = %.10f\n', F.doping);
fprintf('filling   = %.10f\n', F.filling);
fprintf('dx        = %.15e\n', F.dx);
fprintf('dy        = %.15e\n', F.dy);
fprintf('mesh_type = %s\n', F.mesh_type);

% ============================================================
% Convert q range to integer range
% ============================================================

iq_min = round(qx_min / F.dx);
iq_max = round(qx_max / F.dx);
jq_min = round(qy_min / F.dy);
jq_max = round(qy_max / F.dy);

iq_list = iq_min:iq_max;
jq_list = jq_min:jq_max;

Niq = numel(iq_list);
Njq = numel(jq_list);

fprintf('\nq scan integer range:\n');
fprintf('iq = [%d, %d], jq = [%d, %d]\n', iq_min, iq_max, jq_min, jq_max);
fprintf('Nq = %d x %d = %d\n', Niq, Njq, Niq * Njq);

% ============================================================
% Calculate chi map
% ============================================================

chi_map = complex(zeros(Njq, Niq));
nPair_map = zeros(Njq, Niq);

for jj = 1:Njq
    jq_q = jq_list(jj);

    fprintf('jq %d / %d : jq = %d\n', jj, Njq, jq_q);

    for ii = 1:Niq
        iq_q = iq_list(ii);

        [chi, nPair] = cal_chi_single_q_from_bin(F, iq_q, jq_q, eta, ...
            use_form_factor, boundary_periodic);

        chi_map(jj, ii) = chi;
        nPair_map(jj, ii) = nPair;
    end
end

qx_list = iq_list * F.dx;
qy_list = jq_list * F.dy;

% ============================================================
% Find max nesting q
% ============================================================

switch plot_quantity
    case "real"
        Qmap = real(chi_map);
        max_label = "max real(chi)";
    case "imag"
        Qmap = imag(chi_map);
        max_label = "max imag(chi)";
    case "abs"
        Qmap = abs(chi_map);
        max_label = "max abs(chi)";
    otherwise
        error('Unknown plot_quantity');
end

% usually skip q = 0 if included
skip_q0 = true;
Qsearch = Qmap;

if skip_q0
    [~, iq0_idx] = min(abs(iq_list));
    [~, jq0_idx] = min(abs(jq_list));
    if iq_list(iq0_idx) == 0 && jq_list(jq0_idx) == 0
        Qsearch(jq0_idx, iq0_idx) = -inf;
    end
end

[max_val, lin_id] = max(Qsearch(:));
[jmax, imax] = ind2sub(size(Qsearch), lin_id);

iq_best = iq_list(imax);
jq_best = jq_list(jmax);

qx_best = iq_best * F.dx;
qy_best = jq_best * F.dy;

chi_best = chi_map(jmax, imax);

fprintf('\n========================================\n');
fprintf('%s\n', max_label);
fprintf('iq_best = %d\n', iq_best);
fprintf('jq_best = %d\n', jq_best);
fprintf('qx_best = %.15e /b1\n', qx_best);
fprintf('qy_best = %.15e /b2\n', qy_best);
fprintf('chi_real = %.15e\n', real(chi_best));
fprintf('chi_imag = %.15e\n', imag(chi_best));
fprintf('chi_abs  = %.15e\n', abs(chi_best));
fprintf('nPair    = %.0f\n', nPair_map(jmax, imax));
fprintf('area_density = %.15e\n', F.dx * F.dy);
fprintf('========================================\n');

% ============================================================
% Plot chi map
% ============================================================

figure;
imagesc(qx_list, qy_list, Qmap);
set(gca, 'YDir', 'normal');
axis equal tight;
colorbar;
colormap hot;

hold on;
plot(qx_best, qy_best, 'wo', 'MarkerSize', 10, 'LineWidth', 2);
plot(qx_best, qy_best, 'kx', 'MarkerSize', 10, 'LineWidth', 2);

xlabel('q_x in b_1 units');
ylabel('q_y in b_2 units');

title(sprintf('\\chi(q), EF=%.4f eV, best=(%.5f, %.5f)', ...
    F.EF, qx_best, qy_best));

set(gca, 'FontSize', 14, 'LineWidth', 1.2);
box on;

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

function [chi, nPair] = cal_chi_single_q_from_bin(F, iq_q, jq_q, eta, ...
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

area_density = F.dx * F.dy;
chi = chi_sum * area_density;

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