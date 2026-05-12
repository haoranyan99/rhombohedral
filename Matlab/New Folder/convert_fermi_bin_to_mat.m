function convert_fermi_bin_to_mat()
clc; close all;

% =========================================================
% Select binary file
% =========================================================
[fname, fpath] = uigetfile('*.bin', ...
    'Select fermiPatch binary file');

if isequal(fname,0)
    return;
end

binfile = fullfile(fpath, fname);

fprintf('Reading:\n%s\n\n', binfile);

% =========================================================
% Open binary
% =========================================================
fid = fopen(binfile, 'rb');

if fid < 0
    error('Cannot open file.');
end

% =========================================================
% Read header
% =========================================================
magic   = fread(fid, 1, 'int32');
version = fread(fid, 1, 'int32');

NkTot = fread(fid, 1, 'int32');
dim   = fread(fid, 1, 'int32');

EF  = fread(fid, 1, 'double');
T_K = fread(fid, 1, 'double');

fprintf('magic   = %d\n', magic);
fprintf('version = %d\n', version);
fprintf('NkTot   = %d\n', NkTot);
fprintf('dim     = %d\n', dim);
fprintf('EF      = %.8f eV\n', EF);
fprintf('T       = %.4f K\n\n', T_K);

% =========================================================
% Allocate
% =========================================================
iq  = zeros(NkTot,1,'int32');
jq  = zeros(NkTot,1,'int32');

kx = zeros(NkTot,1);
ky = zeros(NkTot,1);

occ_avg = zeros(NkTot,1);

evals    = zeros(NkTot, dim);
occ_band = zeros(NkTot, dim);

% eigenvectors:
% size = [NkTot, dim, dim]
% U(k, band, orbital)
evec = complex( ...
    zeros(NkTot, dim, dim), ...
    zeros(NkTot, dim, dim));

% =========================================================
% Read each k point
% =========================================================
for ik = 1:NkTot

    iq(ik) = fread(fid, 1, 'int32');
    jq(ik) = fread(fid, 1, 'int32');

    kx(ik) = fread(fid, 1, 'double');
    ky(ik) = fread(fid, 1, 'double');

    occ_avg(ik) = fread(fid, 1, 'double');

    evals(ik,:) = fread(fid, dim, 'double');

    occ_band(ik,:) = fread(fid, dim, 'double');

    evec_re = fread(fid, dim*dim, 'double');
    evec_im = fread(fid, dim*dim, 'double');

    tmp = complex(evec_re, evec_im);

    % -----------------------------------------
    % reshape:
    %
    % stored as:
    % evec[b * dim + a]
    %
    % => row = band
    % => col = orbital
    % -----------------------------------------
    tmp = reshape(tmp, [dim, dim]);

    evec(ik,:,:) = tmp;
end

fclose(fid);

fprintf('Finished reading.\n');

% =========================================================
% Build structure
% =========================================================
FERMI = struct();

FERMI.magic   = magic;
FERMI.version = version;

FERMI.NkTot = NkTot;
FERMI.dim   = dim;

FERMI.EF  = EF;
FERMI.T_K = T_K;

FERMI.iq = iq;
FERMI.jq = jq;

FERMI.kx = kx;
FERMI.ky = ky;

FERMI.occ_avg  = occ_avg;
FERMI.evals    = evals;
FERMI.occ_band = occ_band;

FERMI.evec = evec;

% =========================================================
% Save mat
% =========================================================
[outname,~,~] = fileparts(fname);

matfile = fullfile(fpath, [outname '.mat']);

save(matfile, 'FERMI', '-v7.3');

fprintf('\nSaved:\n%s\n', matfile);

end