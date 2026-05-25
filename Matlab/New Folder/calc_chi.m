function calc_chi()

% ============================================================
% USER SETTINGS
% qx,qy are in units of b1,b2:
%   q = qx * b1 + qy * b2
% bin grid:
%   qx = iq * F.dx
%   qy = jq * F.dy
% ============================================================
qx_target = 0.005;   % in /b1
qy_target = -0.005;   % in /b2

eta = 1e-5;
use_form_factor = true;
boundary_periodic = false;

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
% Convert qx,qy -> integer iq,jq
% ============================================================
iq_target = round(qx_target / F.dx);
jq_target = round(qy_target / F.dy);

qx_actual = iq_target * F.dx;
qy_actual = jq_target * F.dy;

% fprintf('\nTarget q fractional:\n');
% fprintf('qx_target = %.15e /b1\n', qx_target);
% fprintf('qy_target = %.15e /b2\n', qy_target);
% 
% fprintf('\nMapped integer q:\n');
% fprintf('iq_target = %d\n', iq_target);
% fprintf('jq_target = %d\n', jq_target);

% ============================================================
% Calculate chi
% ============================================================
chi = cal_chi_single_q_from_bin(F, iq_target, jq_target, eta, ...
    use_form_factor, boundary_periodic);

fprintf('\nq integer = (%d, %d)\n', iq_target, jq_target);
% fprintf('q fractional = (%.15e, %.15e) in b1,b2 units\n', qx_actual, qy_actual);
fprintf('chi_real = %.15e\n', real(chi));
fprintf('chi_imag = %.15e\n', imag(chi));
fprintf('========================================\n');

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

    % C++ flatten:
    % pos = orbital_index * dim + band_index
    %
    % MATLAB reshape is column-major, so transpose is required.
    %
    % After this:
    % evec_re(ik, orbital, band)
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

% ============================================================
% Normalization
% ============================================================
chi_patch_avg = chi_sum;
area_density = F.dx * F.dy;

chi = chi_patch_avg * area_density;

% fprintf('nPair = %d\n', nPair);
% fprintf('effective k pairs = %.6f\n', nPair / (dim * dim));

% fprintf('chi_patch_avg_real = %.15e\n', real(chi_patch_avg));
% fprintf('chi_patch_avg_imag = %.15e\n', imag(chi_patch_avg));
fprintf('area_density       = %.15e\n', area_density);

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