function plot_fermi_contour()

    % --------------------------
    % select file
    % --------------------------
    [fname, fpath] = uigetfile('*.txt', 'Select fermiPatch file');
    if isequal(fname,0); return; end
    file = fullfile(fpath, fname);

    % --------------------------
    % read header: EF as mu + doping
    % --------------------------
    fid = fopen(file, 'r');
    if fid < 0
        error('Cannot open file.');
    end

    mu_val  = NaN;
    dop_val = NaN;

    while ~feof(fid)
        line = fgetl(fid);

        tok = regexp(line, '^\s*#\s*EF\s*=\s*([-\d\.Ee+]+)', 'tokens');
        if ~isempty(tok)
            mu_val = str2double(tok{1}{1});
        end

        tok = regexp(line, '^\s*#\s*doping\s*=\s*([-\d\.Ee+]+)', 'tokens');
        if ~isempty(tok)
            dop_val = str2double(tok{1}{1});
        end

        if ~isnan(mu_val) && ~isnan(dop_val)
            break;
        end
    end

    fclose(fid);

    % --------------------------
    % user params
    % --------------------------
    band_idx = 5;
    thr      = 0.5;

    kx_range = [-1.6, -1.35];
    ky_range = [0.75, 0.95];

    Nx_fine = 500;
    Ny_fine = 500;
    interp_method = 'natural';

    % --------------------------
    % read data
    % --------------------------
    raw = readmatrix(file, "FileType", "text", "CommentStyle", "#");

    iq = raw(:,2);
    jq = raw(:,3);
    kx = raw(:,4);
    ky = raw(:,5);

    dim = (size(raw,2) - 6) / 2;
    fcol = 7 + dim + band_idx;
    f = raw(:,fcol);

    % --------------------------
    % build g1, g2
    % --------------------------
    i00 = find(iq == 0 & jq == 0, 1);
    i10 = find(iq == 1 & jq == 0, 1);
    i01 = find(iq == 0 & jq == 1, 1);

    if isempty(i00) || isempty(i10) || isempty(i01)
        error('Need (iq,jq) = (0,0), (1,0), (0,1) points in file.');
    end

    k00 = [kx(i00); ky(i00)];
    g1  = [kx(i10); ky(i10)] - k00;
    g2  = [kx(i01); ky(i01)] - k00;

    % --------------------------
    % grid from iq,jq
    % --------------------------
    Nk = max(abs([iq; jq]));
    N  = 2 * Nk + 1;

    KX = NaN(N,N);
    KY = NaN(N,N);
    F  = NaN(N,N);

    row = jq + Nk + 1;
    col = iq + Nk + 1;

    for i = 1:length(iq)
        r = row(i);
        c = col(i);

        if isnan(F(r,c))
            kvec = k00 + iq(i) * g1 + jq(i) * g2;
            KX(r,c) = kvec(1);
            KY(r,c) = kvec(2);
            F(r,c)  = f(i);
        end
    end

    % shift origin
    KX = KX - k00(1);
    KY = KY - k00(2);

    % --------------------------
    % crop region
    % --------------------------
    mask = KX >= kx_range(1) & KX <= kx_range(2) & ...
           KY >= ky_range(1) & KY <= ky_range(2);

    valid = isfinite(KX) & isfinite(KY) & isfinite(F) & mask;

    x0 = KX(valid);
    y0 = KY(valid);
    z0 = F(valid);

    % --------------------------
    % interpolate
    % --------------------------
    xfine = linspace(kx_range(1), kx_range(2), Nx_fine);
    yfine = linspace(ky_range(1), ky_range(2), Ny_fine);
    [KXf, KYf] = meshgrid(xfine, yfine);

    Ff = griddata(x0, y0, z0, KXf, KYf, interp_method);

    % --------------------------
    % plot
    % --------------------------
    fig = figure('Color', 'w');
    ax = axes(fig);
    hold(ax, 'on');

    contour(ax, KXf, KYf, Ff, [thr thr], ...
        'k-', 'LineWidth', 1.8);

    xlim(ax, kx_range);
    ylim(ax, ky_range);

    axis(ax, 'equal');
    box(ax, 'on');

    xlabel(ax, 'k_x');
    ylabel(ax, 'k_y');

    title(ax, sprintf('\\mu = %.4f   |   doping = %.4f', ...
        mu_val, dop_val), ...
        'Interpreter', 'tex');

    set(ax, 'FontSize', 14, 'LineWidth', 1.2);

end