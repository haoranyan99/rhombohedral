function make_tiny_fermi_bin_nesting_lines()

clc;

outfile = 'tiny_fermi_nesting_lines.bin';

Nk = 2;
P = 2*Nk + 1;
NkTot = P * P;
dim = 1;

EF = 0.0;
T_K = 1.0;
doping = 0.0;
filling = 0.5;

dx = 1.0 / P;
dy = 1.0 / P;
mesh_type = 'square';

magic = int32(20260510);
version = int32(3);

fid = fopen(outfile, 'wb');
if fid < 0
    error('Cannot open output file.');
end

fwrite(fid, magic, 'int32');
fwrite(fid, version, 'int32');
fwrite(fid, int32(NkTot), 'int32');
fwrite(fid, int32(dim), 'int32');

fwrite(fid, EF, 'double');
fwrite(fid, T_K, 'double');
fwrite(fid, doping, 'double');
fwrite(fid, filling, 'double');

fwrite(fid, dx, 'double');
fwrite(fid, dy, 'double');

mesh_buf = zeros(1, 32, 'uint8');
mesh_bytes = uint8(mesh_type);
mesh_buf(1:numel(mesh_bytes)) = mesh_bytes;
fwrite(fid, mesh_buf, 'uint8');

for iq = -Nk:Nk
    for jq = -Nk:Nk

        kx = double(iq);
        ky = double(jq);

        % Artificial nested two-line model:
        % left line iq=-1: occupied, E=-0.5
        % right line iq=+1: empty, E=+0.5
        % others: empty high energy, E=+2
        if iq == -1
            evals = -0.5;
            occ_band = 1.0;
        elseif iq == 1
            evals = 0.5;
            occ_band = 0.0;
        else
            evals = 2.0;
            occ_band = 0.0;
        end

        occ_avg = occ_band;

        evec_re = 1.0;
        evec_im = 0.0;

        fwrite(fid, int32(iq), 'int32');
        fwrite(fid, int32(jq), 'int32');

        fwrite(fid, kx, 'double');
        fwrite(fid, ky, 'double');
        fwrite(fid, occ_avg, 'double');

        fwrite(fid, evals, 'double');
        fwrite(fid, occ_band, 'double');

        fwrite(fid, evec_re, 'double');
        fwrite(fid, evec_im, 'double');
    end
end

fclose(fid);

fprintf('Generated: %s\n', fullfile(pwd, outfile));
fprintf('NkTot = %d\n', NkTot);
fprintf('dim   = %d\n', dim);
fprintf('\nExpected test:\n');
fprintf('q = (2,0), open boundary, use_form_factor=true\n');
fprintf('|chi_real| = 5/25 = %.15e\n', 5/25);

end