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
    band_idx = 5;        % ← 改这里
    thr      = 0.5;      % ← 改这里（费米面通常 0.5）

    % region (可选裁剪范围)
    kx_range = [-1.6, -1.4];   % [] 表示不裁剪
    ky_range = [0.8, 0.95];

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
    % build g1, g2
    % --------------------------
    i00 = find(iq==0 & jq==0,1);
    i10 = find(iq==1 & jq==0,1);
    i01 = find(iq==0 & jq==1,1);

    k00 = [kx(i00); ky(i00)];
    g1  = [kx(i10); ky(i10)] - k00;
    g2  = [kx(i01); ky(i01)] - k00;

    % --------------------------
    % reconstruct k-space grid
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
    % optional crop region
    % --------------------------
    mask = true(size(KX));

    if ~isempty(kx_range)
        mask = mask & (KX>=kx_range(1) & KX<=kx_range(2));
    end
    if ~isempty(ky_range)
        mask = mask & (KY>=ky_range(1) & KY<=ky_range(2));
    end

    F(~mask) = NaN;

    % --------------------------
    % plot ONLY contour
    % --------------------------
    figure('Color','w'); hold on;

    contour(KX, KY, F, [thr thr], ...
        'k-', 'LineWidth', 2);

    axis equal;
    axis tight;
    box on;

    xlabel('k_x');
    ylabel('k_y');
    title(sprintf('Band %d, f = %.2f contour', band_idx, thr));
    xlim(kx_range);
    ylim(ky_range);

    set(gca,'FontSize',14,'LineWidth',1.2);

end