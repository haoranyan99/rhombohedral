function plot_fermi_scatter_binary()

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
    thr      = 0.5;     % occupation threshold

    % display range: only plot this region
    kx_range = [-1.6, -1.35];
    ky_range = [0.72, 1.00];
%     kx_range = [-2, 2];
%     ky_range = [-2, 2];

    % marker style
    marker_size = 18;

    % colors
    col_empty  = [0, 0, 0];          % black
    col_filled = [0.88, 0.88, 0.88]; % very light gray

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
    % reconstruct k-space using iq,jq topology
    % --------------------------
    KX = k00(1) + iq .* g1(1) + jq .* g2(1);
    KY = k00(2) + iq .* g1(2) + jq .* g2(2);

    % shift origin, same as your old script
    KX = KX - k00(1);
    KY = KY - k00(2);

    % --------------------------
    % crop region
    % --------------------------
    mask = (KX >= kx_range(1) & KX <= kx_range(2) & ...
            KY >= ky_range(1) & KY <= ky_range(2));

    KX = KX(mask);
    KY = KY(mask);
    f  = f(mask);
    iq = iq(mask);
    jq = jq(mask);

    % binary classification
    is_filled = (f >= thr);
    is_empty  = ~is_filled;

    % --------------------------
    % plot
    % --------------------------
    fig = figure('Color','w');
    ax = axes(fig);
    hold(ax, 'on');

    % draw filled first, then black on top
    scatter(ax, KX(is_filled), KY(is_filled), marker_size, ...
        'MarkerFaceColor', col_filled, ...
        'MarkerEdgeColor', col_filled);

    scatter(ax, KX(is_empty), KY(is_empty), marker_size, ...
        'MarkerFaceColor', col_empty, ...
        'MarkerEdgeColor', col_empty);

    xlim(ax, kx_range);
    ylim(ax, ky_range);
    box(ax, 'on');
    axis(ax, 'equal');

    xlabel(ax, 'k_x');
    ylabel(ax, 'k_y');
    title(ax, sprintf('Binary occupation scatter (band %d, thr = %.2f)', band_idx, thr));
    set(ax, 'FontSize', 14, 'LineWidth', 1.2);

    % --------------------------
    % optional click info
    % --------------------------
    S.ax = ax;
    S.KX = KX;
    S.KY = KY;
    S.f  = f;
    S.iq = iq;
    S.jq = jq;
    S.band_idx = band_idx;

    guidata(fig, S);
    set(fig, 'WindowButtonDownFcn', @clickCallback);

    function clickCallback(src, ~)
        S = guidata(src);

        cp = get(S.ax, 'CurrentPoint');
        xq = cp(1,1);
        yq = cp(1,2);

        xl = xlim(S.ax);
        yl = ylim(S.ax);
        if xq < xl(1) || xq > xl(2) || yq < yl(1) || yq > yl(2)
            return;
        end

        dx = S.KX - xq;
        dy = S.KY - yq;
        dist2 = dx.^2 + dy.^2;
        [~, imin] = min(dist2);

        fprintf('iq=%d, jq=%d, kx=%.6f, ky=%.6f, f_%d=%.6f\n', ...
            S.iq(imin), S.jq(imin), ...
            S.KX(imin), S.KY(imin), ...
            S.band_idx, S.f(imin));
    end

end