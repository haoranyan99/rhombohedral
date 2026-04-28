function out = phase_diagram_LatticePhaseTrans()
% phase_diagram_lattice_manual_range
% Lattice-only phase diagram with manual doping/T range.
%
% Boundary:
%   delta = (5/6)*b2^2 - a2*c2 = 0
%
% Click any point to show:
%   doping, T, a2, b2, c2, delta

    % =========================================================
    % 1. user settings
    % =========================================================
    dop_min = -2.0;
    dop_max =  0.0;
    Ndop    = 401;

    T_min = 0.0;
    T_max = 10.0;
    NT    = 401;

    show_background = true;
    fs = 18;

    % =========================================================
    % 2. load params / coeff
    % =========================================================
    par  = make_realchi_params(true);
    coef = make_realchi_coeff(par);

    % =========================================================
    % 3. build grid
    % =========================================================
    dop_list = linspace(dop_min, dop_max, Ndop);
    T_list   = linspace(T_min, T_max, NT);

    [DOP, TT] = meshgrid(dop_list, T_list);

    a2_map    = nan(NT, Ndop);
    b2_map    = nan(NT, Ndop);
    c2_map    = nan(NT, Ndop);
    delta_map = nan(NT, Ndop);
    phase_map = nan(NT, Ndop);   % 1 or 2

    % =========================================================
    % 4. evaluate lattice coefficients
    % =========================================================
    for i = 1:NT
        for j = 1:Ndop
            T   = TT(i,j);
            dop = DOP(i,j);

            a2 = coef.a2(T);
            b2 = coef.b2(dop);
            c2 = coef.c2;

            if ~isfinite(a2) || ~isfinite(b2) || ~isfinite(c2)
                continue;
            end

            delta = (5/6) * b2.^2 - a2 .* c2;

            a2_map(i,j)    = a2;
            b2_map(i,j)    = b2;
            c2_map(i,j)    = c2;
            delta_map(i,j) = delta;

            if delta >= 0
                phase_map(i,j) = 1;
            else
                phase_map(i,j) = 2;
            end
        end
    end

    fprintf('valid points = %d / %d\n', nnz(isfinite(delta_map)), numel(delta_map));

    % =========================================================
    % 5. plot
    % =========================================================
    fig = figure('Color','w', 'Name','Lattice phase diagram');
    ax = axes(fig);
    hold(ax, 'on');

    img_handle = [];

    if show_background
        % phase 1 -> red, phase 2 -> blue
        Cplot = nan(size(phase_map));
        Cplot(phase_map == 1) = 1;
        Cplot(phase_map == 2) = 2;

        img_handle = imagesc(ax, dop_list, T_list, Cplot);
        set(ax, 'YDir', 'normal');

        colormap(ax, [0.85 0.20 0.15;
                      0.20 0.45 0.90]);
        caxis(ax, [1 2]);

        cb = colorbar(ax);
        cb.Ticks = [1.25, 1.75];
        cb.TickLabels = {'\Delta \ge 0', '\Delta < 0'};
        cb.FontSize = fs;
    end

    % boundary delta = 0
    contour_handle = contour(ax, dop_list, T_list, delta_map, [0 0], ...
        'k-', 'LineWidth', 2.0);

    xlabel(ax, 'doping', 'FontSize', fs);
    ylabel(ax, 'T', 'FontSize', fs);
    set(ax, 'FontSize', fs, 'LineWidth', 1.2, 'Box', 'on');
    axis(ax, 'tight');

    % =========================================================
    % 6. data tip
    % =========================================================
    dcm = datacursormode(fig);
    set(dcm, 'Enable', 'on', ...
             'SnapToDataVertex', 'off', ...
             'DisplayStyle', 'datatip');
    set(dcm, 'UpdateFcn', @tip_cb_);

    if isgraphics(img_handle)
        set(img_handle, 'PickableParts', 'all', 'HitTest', 'on');
    end
    if isgraphics(contour_handle)
        set(contour_handle, 'PickableParts', 'all', 'HitTest', 'on');
    end

    % =========================================================
    % 7. output
    % =========================================================
    out = struct();
    out.dop_list   = dop_list;
    out.T_list     = T_list;
    out.a2_map     = a2_map;
    out.b2_map     = b2_map;
    out.c2_map     = c2_map;
    out.delta_map  = delta_map;
    out.phase_map  = phase_map;

    % =========================================================
    % nested datatip callback
    % =========================================================
    function txt = tip_cb_(~, evt)
        pos = evt.Position;
        x = pos(1);   % doping
        y = pos(2);   % T

        [~, j] = min(abs(dop_list - x));
        [~, i] = min(abs(T_list   - y));

        dop0 = dop_list(j);
        T0   = T_list(i);

        a20 = a2_map(i,j);
        b20 = b2_map(i,j);
        c20 = c2_map(i,j);
        d0  = delta_map(i,j);
        ph  = phase_map(i,j);

        if isnan(d0)
            txt = {
                sprintf('doping = %.6g', dop0)
                sprintf('T = %.6g', T0)
                'invalid point'
            };
            return;
        end

        txt = {
            sprintf('doping = %.6g', dop0)
            sprintf('T = %.6g', T0)
            sprintf('phase = %d', ph)
            sprintf('a2 = %.6g', a20)
            sprintf('b2 = %.6g', b20)
            sprintf('c2 = %.6g', c20)
            sprintf('delta = %.6g', d0)
            'delta = 5/6*b2^2 - a2*c2'
        };
    end
end