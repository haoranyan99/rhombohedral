function plot_fermi_contour()

    clc; close all;

    % --------------------------
    % select file
    % --------------------------
    [fname, fpath] = uigetfile('*.txt', 'Select fermiPatch file');
    if isequal(fname,0); return; end
    file = fullfile(fpath, fname);

    % --------------------------
    % user params
    % --------------------------
    band_idx = 5;
    thr      = 0.5;

    % display range
    % kx_range = [-1.55, -1.4];
    % ky_range = [0.8, 0.92];
    kx_range = [-1.6, -1.35];
    ky_range = [0.72, 1];

    % smoothing grid
    Nx_fine = 500;
    Ny_fine = 500;
    interp_method = 'natural';   % 'nearest' / 'linear' / 'natural' / 'cubic' / 'v4'

    % --------------------------
    % read data
    % --------------------------
    raw = readmatrix(file, "FileType","text", "CommentStyle","#");

    iq = raw(:,2);
    jq = raw(:,3);
    kx = raw(:,4);
    ky = raw(:,5);

    dim = (size(raw,2)-6)/2;
    fcol = 7 + dim + band_idx;
    f = raw(:,fcol);

    % --------------------------
    % build g1, g2 from (0,0),(1,0),(0,1)
    % --------------------------
    i00 = find(iq==0 & jq==0,1);
    i10 = find(iq==1 & jq==0,1);
    i01 = find(iq==0 & jq==1,1);

    if isempty(i00) || isempty(i10) || isempty(i01)
        error('Need (iq,jq)=(0,0),(1,0),(0,1) points in file.');
    end

    k00 = [kx(i00); ky(i00)];
    g1  = [kx(i10); ky(i10)] - k00;
    g2  = [kx(i01); ky(i01)] - k00;

    % --------------------------
    % grid from iq,jq topology
    % --------------------------
    Nk = max(abs([iq; jq]));
    N  = 2*Nk + 1;

    KX = NaN(N,N);
    KY = NaN(N,N);
    F  = NaN(N,N);

    row = jq + Nk + 1;
    col = iq + Nk + 1;

    for i = 1:length(iq)
        r = row(i); c = col(i);
        if isnan(F(r,c))
            kvec = k00 + iq(i)*g1 + jq(i)*g2;
            KX(r,c) = kvec(1);
            KY(r,c) = kvec(2);
            F(r,c)  = f(i);
        end
    end

    % shift origin
    KX = KX - k00(1);
    KY = KY - k00(2);

    % --------------------------
    % raw points for nearest-point query
    % --------------------------
    KX_pts = KX(:);
    KY_pts = KY(:);
    F_pts  = F(:);

    IQ_pts = repmat((-Nk:Nk), N, 1);
    IQ_pts = IQ_pts(:);

    JQ_pts = repmat((-Nk:Nk)', 1, N);
    JQ_pts = JQ_pts(:);

    valid = isfinite(KX_pts) & isfinite(KY_pts) & isfinite(F_pts);
    KX_pts = KX_pts(valid);
    KY_pts = KY_pts(valid);
    F_pts  = F_pts(valid);
    IQ_pts = IQ_pts(valid);
    JQ_pts = JQ_pts(valid);

    in_region = (KX_pts>=kx_range(1) & KX_pts<=kx_range(2) & ...
                 KY_pts>=ky_range(1) & KY_pts<=ky_range(2));

    KX_pts = KX_pts(in_region);
    KY_pts = KY_pts(in_region);
    F_pts  = F_pts(in_region);
    IQ_pts = IQ_pts(in_region);
    JQ_pts = JQ_pts(in_region);

    % --------------------------
    % crop original grid first
    % --------------------------
    mask = (KX>=kx_range(1) & KX<=kx_range(2) & ...
            KY>=ky_range(1) & KY<=ky_range(2));

    KXc = KX;
    KYc = KY;
    Fc  = F;

    KXc(~mask) = NaN;
    KYc(~mask) = NaN;
    Fc(~mask)  = NaN;

    % --------------------------
    % smooth contour by interpolation onto fine grid
    % --------------------------
    valid2 = isfinite(KXc) & isfinite(KYc) & isfinite(Fc);
    x0 = KXc(valid2);
    y0 = KYc(valid2);
    z0 = Fc(valid2);

    xfine = linspace(kx_range(1), kx_range(2), Nx_fine);
    yfine = linspace(ky_range(1), ky_range(2), Ny_fine);
    [KXf, KYf] = meshgrid(xfine, yfine);

    Ff = griddata(x0, y0, z0, KXf, KYf, interp_method);

    % --------------------------
    % plot
    % --------------------------
    fig = figure('Color','w');
    ax = axes(fig);
    hold(ax, 'on');

    contour(ax, KXf, KYf, Ff, [thr thr], 'k-', 'LineWidth', 1.8);

    xlim(ax, kx_range);
    ylim(ax, ky_range);
    box(ax, 'on');

    xlabel(ax, 'k_x');
    ylabel(ax, 'k_y');
    title(ax, sprintf('Band %d, f = %.2f contour', band_idx, thr));
    set(ax, 'FontSize', 14, 'LineWidth', 1.2);

    % --------------------------
    % store data
    % --------------------------
    S.ax = ax;
    S.KX_pts = KX_pts;
    S.KY_pts = KY_pts;
    S.F_pts  = F_pts;
    S.IQ_pts = IQ_pts;
    S.JQ_pts = JQ_pts;
    S.band_idx = band_idx;
    S.marker = [];
    S.label  = [];
    guidata(fig, S);

    % left click = show/update label
    % right click = clear label
    set(fig, 'WindowButtonDownFcn', @mouseClickCallback);

    function mouseClickCallback(src, ~)
        S = guidata(src);

        % selection type:
        % normal  = left click
        % alt     = right click / ctrl-click
        clickType = get(src, 'SelectionType');

        if strcmp(clickType, 'alt')
            if isgraphics(S.marker), delete(S.marker); end
            if isgraphics(S.label),  delete(S.label);  end
            S.marker = [];
            S.label  = [];
            guidata(src, S);
            return;
        end

        cp = get(S.ax, 'CurrentPoint');
        xq = cp(1,1);
        yq = cp(1,2);

        xl = xlim(S.ax);
        yl = ylim(S.ax);
        if xq < xl(1) || xq > xl(2) || yq < yl(1) || yq > yl(2)
            return;
        end

        dx = S.KX_pts - xq;
        dy = S.KY_pts - yq;
        dist2 = dx.^2 + dy.^2;
        [~, imin] = min(dist2);

        iq_near = S.IQ_pts(imin);
        jq_near = S.JQ_pts(imin);
        kx_near = S.KX_pts(imin);
        ky_near = S.KY_pts(imin);
        f_near  = S.F_pts(imin);

        if isgraphics(S.marker), delete(S.marker); end
        if isgraphics(S.label),  delete(S.label);  end

        S.marker = plot(S.ax, kx_near, ky_near, 'ro', ...
            'MarkerSize', 7, 'LineWidth', 1.3);

        S.label = text(S.ax, kx_near, ky_near, ...
            sprintf('  (iq,jq)=(%d,%d)\n  f_%d=%.4f', ...
            iq_near, jq_near, S.band_idx, f_near), ...
            'Color', 'r', 'FontSize', 11, ...
            'VerticalAlignment', 'bottom', ...
            'Interpreter', 'none');

        guidata(src, S);

        fprintf('nearest point: iq=%d, jq=%d, kx=%.6f, ky=%.6f, f_%d=%.6f\n', ...
            iq_near, jq_near, kx_near, ky_near, S.band_idx, f_near);
    end
end