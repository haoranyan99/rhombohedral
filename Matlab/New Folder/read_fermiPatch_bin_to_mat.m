function F = plot_fermiPatch_bin_energy_clicktip()
% ============================================================
% Read fermiPatch_*.bin and plot energy scatter with custom click tips
%
% 手写参数在这里改：
%   band_idx        : 要画的 band，1~dim
%   plot_E_minus_EF : true 画 E-EF；false 画 E
%   energy_unit     : 'eV' or 'meV'
%   marker_size     : 散点大小
%   near_EF_meV     : 是否额外标记 EF 附近点，设 [] 关闭
%
% 点击一个点：显示 Energy 和 occ
% 再次点击同一个点：取消该点提示
% 可同时显示多个点提示
% ============================================================

    clearvars -except F;

    % ========================================================
    % USER SETTINGS
    % ========================================================
    band_idx = 6;              % 手动改这里：1 ~ dim
    plot_E_minus_EF = false;   % true: color = E-EF, false: color = E
    energy_unit = 'eV';        % 'eV' or 'meV'

    marker_size = 22;
    marker_alpha = 0.9;

    colormap_name = 'turbo';

    near_EF_meV = [];          % 例如 2 表示标记 |E-EF|<2 meV；[] 表示不标记

    title_fontsize = 14;
    axis_fontsize  = 13;

    % ========================================================
    % Select bin file only
    % ========================================================
    [fname, fpath] = uigetfile('*.bin', 'Select fermiPatch bin file');

    if isequal(fname, 0)
        disp('Canceled.');
        F = [];
        return;
    end

    bin_file = fullfile(fpath, fname);

    % ========================================================
    % Read bin
    % ========================================================
    F = read_fermiPatch_bin_no_evec(bin_file);
    assignin('base', 'F', F);

    if band_idx < 1 || band_idx > F.dim
        error('band_idx = %d is invalid. This file has dim = %d bands.', band_idx, F.dim);
    end

    % ========================================================
    % Prepare data
    % ========================================================
    Eraw = F.evals(:, band_idx);
    occ  = F.occ_band(:, band_idx);

    switch lower(energy_unit)
        case 'ev'
            unit_factor = 1;
            unit_label = 'eV';
        case 'mev'
            unit_factor = 1000;
            unit_label = 'meV';
        otherwise
            error('energy_unit must be ''eV'' or ''meV''.');
    end

    if plot_E_minus_EF
        C = unit_factor * (Eraw - F.EF);
        color_label = sprintf('E - E_F (%s)', unit_label);
    else
        C = unit_factor * Eraw;
        color_label = sprintf('Energy (%s)', unit_label);
    end

    % ========================================================
    % Print summary
    % ========================================================
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('File      : %s\n', bin_file);
    fprintf('NkTot     : %d\n', F.NkTot);
    fprintf('dim       : %d\n', F.dim);
    fprintf('band_idx  : %d\n', band_idx);
    fprintf('EF        : %.12f eV\n', F.EF);
    fprintf('T_K       : %.6f K\n', F.T_K);
    fprintf('kBT       : %.6f meV\n', 8.617333262e-2 * F.T_K);
    fprintf('doping    : %.12f\n', F.doping);
    fprintf('mesh_type : %s\n', F.mesh_type);
    fprintf('E range   : [%.8f, %.8f] eV\n', min(Eraw), max(Eraw));
    fprintf('occ range : [%.8f, %.8f]\n', min(occ), max(occ));
    fprintf('============================================================\n');

    % ========================================================
    % Plot
    % ========================================================
    fig = figure('Color', 'w');
    ax = axes(fig);
    hold(ax, 'on');

    s = scatter(ax, F.kx, F.ky, marker_size, C, 'filled');
    s.MarkerFaceAlpha = marker_alpha;
    s.MarkerEdgeAlpha = marker_alpha;

    axis(ax, 'equal');
    box(ax, 'on');
    xlabel(ax, 'k_x (1/A)');
    ylabel(ax, 'k_y (1/A)');

    colormap(ax, colormap_name);
    cb = colorbar(ax);
    ylabel(cb, color_label);

    title(ax, sprintf('Band %d energy distribution, EF = %.6f eV', band_idx, F.EF), ...
        'FontSize', title_fontsize);

    set(ax, 'FontSize', axis_fontsize);
    set(ax, 'LineWidth', 1.2);
    set(ax, 'TickDir', 'in');

    % Optional EF marker
    if ~isempty(near_EF_meV)
        idx_near = abs(1000 * (Eraw - F.EF)) < near_EF_meV;
        if any(idx_near)
            plot(ax, F.kx(idx_near), F.ky(idx_near), 'k.', 'MarkerSize', 6);
        end
    end

    % ========================================================
    % Custom click tips
    % ========================================================
    tip_data.kx = F.kx;
    tip_data.ky = F.ky;
    tip_data.Eraw = Eraw;
    tip_data.occ = occ;
    tip_data.EF = F.EF;
    tip_data.band_idx = band_idx;
    tip_data.unit_factor = unit_factor;
    tip_data.unit_label = unit_label;
    tip_data.active_idx = [];
    tip_data.text_handles = gobjects(0);

    setappdata(fig, 'tip_data', tip_data);

    s.ButtonDownFcn = @(src, event) toggle_point_tip(src, event, fig, ax);

    fprintf('\nClick a point to show Energy and occ.\n');
    fprintf('Click the same point again to remove its label.\n');
    fprintf('Multiple labels can stay on the figure.\n');

end


% ============================================================
% Read bin, skip eigenvectors
% ============================================================
function F = read_fermiPatch_bin_no_evec(bin_file)

    fid = fopen(bin_file, 'rb');
    if fid < 0
        error('Cannot open file: %s', bin_file);
    end
    cleaner = onCleanup(@() fclose(fid));

    F.magic   = fread(fid, 1, 'int32');
    F.version = fread(fid, 1, 'int32');

    F.NkTot   = fread(fid, 1, 'int32');
    F.dim     = fread(fid, 1, 'int32');

    F.EF      = fread(fid, 1, 'double');
    F.T_K     = fread(fid, 1, 'double');
    F.doping  = fread(fid, 1, 'double');
    F.filling = fread(fid, 1, 'double');

    F.dx      = fread(fid, 1, 'double');
    F.dy      = fread(fid, 1, 'double');

    mesh_raw = fread(fid, 32, '*char')';
    F.mesh_type = string(regexprep(mesh_raw, char(0), ''));

    if F.magic ~= 20260510
        error('Wrong magic number: %d. This may not be fermiPatch bin.', F.magic);
    end

    if F.version ~= 3
        warning('Expected version 3, but got version %d.', F.version);
    end

    NkTot = double(F.NkTot);
    dim   = double(F.dim);

    F.iq      = zeros(NkTot, 1);
    F.jq      = zeros(NkTot, 1);
    F.kx      = zeros(NkTot, 1);
    F.ky      = zeros(NkTot, 1);
    F.occ_avg = zeros(NkTot, 1);

    F.evals    = zeros(NkTot, dim);
    F.occ_band = zeros(NkTot, dim);

    for ik = 1:NkTot

        F.iq(ik) = fread(fid, 1, 'int32');
        F.jq(ik) = fread(fid, 1, 'int32');

        F.kx(ik) = fread(fid, 1, 'double');
        F.ky(ik) = fread(fid, 1, 'double');

        F.occ_avg(ik) = fread(fid, 1, 'double');

        tmpE = fread(fid, dim, 'double');
        tmpF = fread(fid, dim, 'double');

        if numel(tmpE) ~= dim || numel(tmpF) ~= dim
            error('Unexpected EOF while reading evals/occ at ik = %d.', ik);
        end

        F.evals(ik, :)    = tmpE.';
        F.occ_band(ik, :) = tmpF.';

        % skip evec_re and evec_im
        offset_bytes = 2 * dim * dim * 8;
        status = fseek(fid, offset_bytes, 'cof');

        if status ~= 0
            error('Unexpected EOF while skipping eigenvectors at ik = %d.', ik);
        end
    end

end


% ============================================================
% Toggle click tip
% ============================================================
function toggle_point_tip(~, event, fig, ax)

    tip_data = getappdata(fig, 'tip_data');

    click_pos = event.IntersectionPoint;
    x0 = click_pos(1);
    y0 = click_pos(2);

    x = tip_data.kx;
    y = tip_data.ky;

    % distance in data coordinate
    d2 = (x - x0).^2 + (y - y0).^2;
    [~, idx] = min(d2);

    % If this point already has a label, delete it
    hit = find(tip_data.active_idx == idx, 1);

    if ~isempty(hit)
        if isgraphics(tip_data.text_handles(hit))
            delete(tip_data.text_handles(hit));
        end
        tip_data.active_idx(hit) = [];
        tip_data.text_handles(hit) = [];
        setappdata(fig, 'tip_data', tip_data);
        return;
    end

    % Otherwise create a new label
    E_eV = tip_data.Eraw(idx);
    occ  = tip_data.occ(idx);

    E_display = tip_data.unit_factor * E_eV;

    label_str = sprintf('E = %.8f %s\nocc = %.8f', ...
        E_display, tip_data.unit_label, occ);

    % small offset so text does not cover the point exactly
    xl = xlim(ax);
    yl = ylim(ax);
    dx = 0.015 * (xl(2) - xl(1));
    dy = 0.015 * (yl(2) - yl(1));

    htxt = text(ax, x(idx) + dx, y(idx) + dy, label_str, ...
        'FontSize', 11, ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k', ...
        'Margin', 4, ...
        'Interpreter', 'none', ...
        'Clipping', 'on');

    tip_data.active_idx(end+1) = idx;
    tip_data.text_handles(end+1) = htxt;

    setappdata(fig, 'tip_data', tip_data);

end