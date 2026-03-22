function out = phase_diagram_four_phases_el_lat()
% phase_diagram_four_phases_el_lat
%
% 4-phase diagram from:
%   electronic order: on/off
%   lattice order:    on/off
%
% Electronic rule (following current electronic script):
%   a1 > 0  -> electronic OFF
%   a1 < 0  -> electronic ON
%
% Lattice rule (following current lattice script):
%   delta = (5/6)*b2^2 - a2*c2
%   delta >= 0 -> lattice ON
%   delta <  0 -> lattice OFF
%
% 4 phases:
%   1 = el off, lat off
%   2 = el on,  lat off
%   3 = el off, lat on
%   4 = el on,  lat on

    % =========================================================
    % 1. params / coeff / style
    % =========================================================
    par  = make_realchi_params(true);
    coef = make_realchi_coeff(par);

    fs = 18;

    % rendered doping grid
    dop_min_user = -4.0;
    dop_max_user =  0.0;
    Ndop_render  = 401;

    % nearest-neighbor tolerance for electronic rendering
    tol_factor = 0.1;

    % whether to overlay boundaries
    show_el_boundary  = true;   % a1 = 0
    show_lat_boundary = true;   % delta = 0

    % =========================================================
    % 2. choose folder
    % =========================================================
    default_root = "/Users/haoranyan/rg_master/data/";
    if ~isfolder(default_root)
        default_root = string(pwd);
    end

    root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
    if isequal(root, 0)
        error('User cancelled.');
    end
    root = string(root);

    iq_pick = round(par.iq_pick);
    jq_pick = round(par.jq_pick);

    fprintf("Root folder:\n  %s\n", root);
    fprintf("Pick absolute q-point: (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);

    % =========================================================
    % 3. load scattered electronic data
    % =========================================================
    dbg = struct();
    dbg.enable = true;
    dbg.max_print = 60;

    [S, dbg_out] = load_chi_scattered_T_doping_debug_(root, iq_pick, jq_pick, dbg);

    if dbg.enable
        print_debug_summary_scattered_(dbg_out);
    end

    if isempty(S.T) || isempty(S.doping) || isempty(S.chi)
        error('No valid scattered (T, doping, chi) points were found.');
    end

    % =========================================================
    % 4. classify raw scattered points by a1 sign
    % =========================================================
    Nraw = numel(S.T);
    a1_raw     = nan(Nraw,1);
    el_on_raw  = nan(Nraw,1);   % 0=off, 1=on

    for k = 1:Nraw
        T   = S.T(k);
        dop = S.doping(k);
        chi_used = S.chi(k);

        try
            C = coef.eval(T, dop, chi_used);
            a1 = C.a1;
        catch
            continue;
        end

        if ~isfinite(a1)
            continue;
        end

        a1_raw(k) = a1;

        % electronic OFF if a1>0 ; ON if a1<0
        if a1 > 0
            el_on_raw(k) = 0;
        elseif a1 < 0
            el_on_raw(k) = 1;
        else
            el_on_raw(k) = NaN;
        end
    end

    valid_mask = isfinite(S.T) & isfinite(S.doping) & isfinite(a1_raw) & isfinite(el_on_raw);

    T_raw      = S.T(valid_mask);
    dop_raw    = S.doping(valid_mask);
    chi_raw    = S.chi(valid_mask);
    a1_raw     = a1_raw(valid_mask);
    el_on_raw  = el_on_raw(valid_mask);
    file_raw   = S.file(valid_mask);

    if isempty(T_raw)
        error('All raw points became invalid after evaluating a1.');
    end

    % =========================================================
    % 5. build common rendered (T, doping) grid from electronic data
    %    first render electronic phase continuously row-by-row
    % =========================================================
    T_list = sort(unique(T_raw));

    if isempty(dop_min_user) || isempty(dop_max_user) || ~isfinite(dop_min_user) || ~isfinite(dop_max_user)
        dop_min = min(dop_raw);
        dop_max = max(dop_raw);
    else
        dop_min = dop_min_user;
        dop_max = dop_max_user;
    end

    dop_list = linspace(dop_min, dop_max, Ndop_render);

    NT   = numel(T_list);
    Ndop = numel(dop_list);

    % rendered electronic maps
    el_on_map   = nan(NT, Ndop);   % 0/1
    a1_map      = nan(NT, Ndop);
    chi_map     = nan(NT, Ndop);
    dop_src_map = nan(NT, Ndop);
    dist_map    = nan(NT, Ndop);
    file_map    = strings(NT, Ndop);

    % do not overfill very large empty gaps
    max_gap_factor = 3.0;

    for iT = 1:NT
        T0 = T_list(iT);
        mT = (T_raw == T0);

        dop_row   = dop_raw(mT);
        a1_row    = a1_raw(mT);
        chi_row   = chi_raw(mT);
        el_row    = el_on_raw(mT);
        file_row  = file_raw(mT);

        if isempty(dop_row)
            continue;
        end

        % sort by doping
        [dop_row, order] = sort(dop_row);
        a1_row   = a1_row(order);
        chi_row  = chi_row(order);
        el_row   = el_row(order);
        file_row = file_row(order);

        % merge duplicated doping values
        [dop_u, ~, ic] = unique(round(dop_row, 10));
        a1_u   = nan(size(dop_u));
        chi_u  = nan(size(dop_u));
        el_u   = nan(size(dop_u));
        file_u = strings(size(dop_u));

        for k = 1:numel(dop_u)
            mm = (ic == k);

            a1_u(k)  = mean(a1_row(mm), 'omitnan');
            chi_u(k) = mean(chi_row(mm), 'omitnan');

            vals = el_row(mm);
            vals = vals(isfinite(vals));
            if ~isempty(vals)
                el_u(k) = round(mean(vals));   % 0 or 1
            end

            idx_first = find(mm, 1, 'first');
            if ~isempty(idx_first)
                file_u(k) = file_row(idx_first);
            end
        end

        dop_u = double(dop_u);

        % keep only good points
        good = isfinite(dop_u) & isfinite(a1_u) & isfinite(chi_u) & isfinite(el_u);
        dop_u  = dop_u(good);
        a1_u   = a1_u(good);
        chi_u  = chi_u(good);
        el_u   = el_u(good);
        file_u = file_u(good);

        if isempty(dop_u)
            continue;
        end

        if numel(dop_u) == 1
            [~, j0] = min(abs(dop_list - dop_u));
            el_on_map(iT, j0)   = el_u;
            a1_map(iT, j0)      = a1_u;
            chi_map(iT, j0)     = chi_u;
            dop_src_map(iT, j0) = dop_u;
            dist_map(iT, j0)    = abs(dop_list(j0) - dop_u);
            file_map(iT, j0)    = file_u;
            continue;
        end

        % midpoint bins
        edges = zeros(1, numel(dop_u)+1);
        edges(2:end-1) = 0.5 * (dop_u(1:end-1) + dop_u(2:end));

        left_gap  = dop_u(2)   - dop_u(1);
        right_gap = dop_u(end) - dop_u(end-1);
        edges(1)   = dop_u(1)   - 0.5 * left_gap;
        edges(end) = dop_u(end) + 0.5 * right_gap;

        raw_spacing_med = median(diff(dop_u));
        render_spacing  = mean(diff(dop_list));

        edges(1)   = max(edges(1), dop_list(1));
        edges(end) = min(edges(end), dop_list(end));

        for k = 1:numel(dop_u)
            bin_width = edges(k+1) - edges(k);

            if bin_width > max_gap_factor * max(raw_spacing_med, render_spacing)
                continue;
            end

            if k < numel(dop_u)
                mask_fill = (dop_list >= edges(k)) & (dop_list < edges(k+1));
            else
                mask_fill = (dop_list >= edges(k)) & (dop_list <= edges(k+1));
            end

            if ~any(mask_fill)
                continue;
            end

            el_on_map(iT, mask_fill)   = el_u(k);
            a1_map(iT, mask_fill)      = a1_u(k);
            chi_map(iT, mask_fill)     = chi_u(k);
            dop_src_map(iT, mask_fill) = dop_u(k);
            dist_map(iT, mask_fill)    = abs(dop_list(mask_fill) - dop_u(k));
            file_map(iT, mask_fill)    = file_u(k);
        end
    end

    % fill single-pixel horizontal holes
    for iT = 1:NT
        row = el_on_map(iT,:);
        for j = 2:(Ndop-1)
            if isnan(row(j)) && isfinite(row(j-1)) && isfinite(row(j+1)) && row(j-1)==row(j+1)
                row(j) = row(j-1);

                if isfinite(a1_map(iT,j-1)) && isfinite(a1_map(iT,j+1))
                    a1_map(iT,j) = 0.5 * (a1_map(iT,j-1) + a1_map(iT,j+1));
                end
                if isfinite(chi_map(iT,j-1)) && isfinite(chi_map(iT,j+1))
                    chi_map(iT,j) = 0.5 * (chi_map(iT,j-1) + chi_map(iT,j+1));
                end
                if isfinite(dop_src_map(iT,j-1)) && isfinite(dop_src_map(iT,j+1))
                    dop_src_map(iT,j) = 0.5 * (dop_src_map(iT,j-1) + dop_src_map(iT,j+1));
                end
                if isfinite(dist_map(iT,j-1)) && isfinite(dist_map(iT,j+1))
                    dist_map(iT,j) = min(dist_map(iT,j-1), dist_map(iT,j+1));
                end

                if file_map(iT,j-1) ~= ""
                    file_map(iT,j) = file_map(iT,j-1);
                else
                    file_map(iT,j) = file_map(iT,j+1);
                end
            end
        end
        el_on_map(iT,:) = row;
    end

    % =========================================================
    % 6. evaluate lattice phase on the SAME rendered grid
    % =========================================================
    a2_map     = nan(NT, Ndop);
    b2_map     = nan(NT, Ndop);
    c2_map     = nan(NT, Ndop);
    delta_map  = nan(NT, Ndop);
    lat_on_map = nan(NT, Ndop);   % 0=off, 1=on

    for i = 1:NT
        for j = 1:Ndop
            T   = T_list(i);
            dop = dop_list(j);

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

            % lattice ON if delta >= 0 ; OFF if delta < 0
            if delta >= 0
                lat_on_map(i,j) = 1;
            else
                lat_on_map(i,j) = 0;
            end
        end
    end

    % =========================================================
    % 7. combine into 4 phases
    % =========================================================
    phase4_map = nan(NT, Ndop);

    for i = 1:NT
        for j = 1:Ndop
            e = el_on_map(i,j);
            l = lat_on_map(i,j);

            if ~isfinite(e) || ~isfinite(l)
                continue;
            end

            if e==0 && l==0
                phase4_map(i,j) = 1;   % el off, lat off
            elseif e==1 && l==0
                phase4_map(i,j) = 2;   % el on, lat off
            elseif e==0 && l==1
                phase4_map(i,j) = 3;   % el off, lat on
            elseif e==1 && l==1
                phase4_map(i,j) = 4;   % el on, lat on
            end
        end
    end

    fprintf('Rendered valid 4-phase points = %d / %d\n', ...
        nnz(isfinite(phase4_map)), numel(phase4_map));

    % =========================================================
    % 8. plot
    % =========================================================
    fig = figure('Color','w', 'Name','Four-phase diagram: electronic + lattice');
    ax = axes(fig);
    hold(ax, 'on');

    Cplot = phase4_map;
    img_handle = imagesc(ax, dop_list, T_list, Cplot);
    set(ax, 'YDir', 'normal');

    % 4 colors
    % 1 = el off, lat off
    % 2 = el on,  lat off
    % 3 = el off, lat on
    % 4 = el on,  lat on
    colormap(ax, [ ...
        0.85 0.85 0.85;   % gray
        0.20 0.45 0.90;   % blue
        0.90 0.55 0.15;   % orange
        0.70 0.20 0.75]); % purple
    caxis(ax, [1 4]);

    cb = colorbar(ax);
    cb.Ticks = [1.375, 2.125, 2.875, 3.625];
    cb.TickLabels = { ...
        'el off, lat off', ...
        'el on, lat off', ...
        'el off, lat on', ...
        'el on, lat on'};
    cb.FontSize = fs;

    if show_el_boundary && any(isfinite(a1_map(:)))
        contour(ax, dop_list, T_list, a1_map, [0 0], 'k--', 'LineWidth', 1.6);
    end

    if show_lat_boundary && any(isfinite(delta_map(:)))
        contour(ax, dop_list, T_list, delta_map, [0 0], 'k-', 'LineWidth', 2.0);
    end

    xlabel(ax, 'doping (10^{12} cm^{-2})', 'FontSize', fs);
    ylabel(ax, 'T (K)', 'FontSize', fs);
    title(ax, '4 phases: electronic on/off + lattice on/off', ...
        'FontSize', fs, 'FontWeight', 'normal');

    set(ax, 'FontSize', fs, 'LineWidth', 1.2, 'Box', 'on');
    axis(ax, 'tight');

    % =========================================================
    % 9. datatip
    % =========================================================
    dcm = datacursormode(fig);
    set(dcm, 'Enable', 'on', ...
             'SnapToDataVertex', 'off', ...
             'DisplayStyle', 'datatip');
    set(dcm, 'UpdateFcn', @tip_cb_);

    if isgraphics(img_handle)
        set(img_handle, 'PickableParts', 'all', 'HitTest', 'on');
    end

    % =========================================================
    % 10. output
    % =========================================================
    out = struct();
    out.dop_list    = dop_list;
    out.T_list      = T_list;
    out.phase4_map  = phase4_map;
    out.el_on_map   = el_on_map;
    out.lat_on_map  = lat_on_map;
    out.a1_map      = a1_map;
    out.a2_map      = a2_map;
    out.b2_map      = b2_map;
    out.c2_map      = c2_map;
    out.delta_map   = delta_map;
    out.chi_map     = chi_map;
    out.dop_src_map = dop_src_map;
    out.dist_map    = dist_map;
    out.file_map    = file_map;
    out.iq_pick     = iq_pick;
    out.jq_pick     = jq_pick;

    % =========================================================
    % nested datatip callback
    % =========================================================
    function txt = tip_cb_(~, evt)
        pos = evt.Position;
        x = pos(1);
        y = pos(2);

        [~, j] = min(abs(dop_list - x));
        [~, i] = min(abs(T_list   - y));

        dop0  = dop_list(j);
        T0    = T_list(i);
        ph    = phase4_map(i,j);
        eon   = el_on_map(i,j);
        lon   = lat_on_map(i,j);
        a10   = a1_map(i,j);
        a20   = a2_map(i,j);
        b20   = b2_map(i,j);
        c20   = c2_map(i,j);
        d0    = delta_map(i,j);
        chi0  = chi_map(i,j);
        dsrc  = dop_src_map(i,j);
        ddist = dist_map(i,j);
        fsrc  = file_map(i,j);

        if isnan(ph)
            txt = {
                sprintf('doping(render) = %.6g', dop0)
                sprintf('T = %.6g', T0)
                'invalid point'
            };
            return;
        end

        txt = {
            sprintf('doping(render) = %.6g', dop0)
            sprintf('T = %.6g', T0)
            sprintf('phase4 = %d', ph)
            sprintf('electronic on = %d', eon)
            sprintf('lattice on = %d', lon)
            sprintf('a1 = %.6g', a10)
            sprintf('chi = %.6g', chi0)
            sprintf('a2 = %.6g', a20)
            sprintf('b2 = %.6g', b20)
            sprintf('c2 = %.6g', c20)
            sprintf('delta = %.6g', d0)
            sprintf('nearest raw doping = %.6g', dsrc)
            sprintf('distance to raw point = %.6g', ddist)
            sprintf('file = %s', fsrc)
        };
    end
end

% =====================================================================
% scattered loader
% =====================================================================
function [S, dbg] = load_chi_scattered_T_doping_debug_(root_dir, iq_pick, jq_pick, dbg_cfg)
    root_dir = string(root_dir);
    L = dir(fullfile(root_dir, "**", "chi*.txt"));
    if isempty(L)
        error("No chi*.txt found under: %s", root_dir);
    end

    iq_pick = round(iq_pick);
    jq_pick = round(jq_pick);

    n = numel(L);

    Tvals = nan(n,1);
    doping_header = nan(n,1);
    chiV = nan(n,1);
    file = strings(n,1);

    n_ok=0; n_fail_header=0; n_fail_read=0; n_fail_cols=0; n_fail_match=0;
    bad_header=strings(0); bad_read=strings(0); bad_cols=strings(0); bad_match=strings(0);

    for k = 1:n
        fpath = string(fullfile(L(k).folder, L(k).name));

        H = parse_header_T_and_doping_robust_(fpath);
        if ~H.ok
            n_fail_header = n_fail_header + 1;
            if dbg_cfg.enable && numel(bad_header) < dbg_cfg.max_print
                bad_header(end+1) = sprintf("%s | header missing T/doping", fpath);
            end
            continue;
        end

        try
            M = read_chi_numeric_skiphash_(fpath);
        catch ME
            n_fail_read = n_fail_read + 1;
            if dbg_cfg.enable && numel(bad_read) < dbg_cfg.max_print
                bad_read(end+1) = sprintf("%s | read failed: %s", fpath, string(ME.message));
            end
            continue;
        end

        if isempty(M) || size(M,2) < 5
            n_fail_cols = n_fail_cols + 1;
            if dbg_cfg.enable && numel(bad_cols) < dbg_cfg.max_print
                bad_cols(end+1) = sprintf("%s | numeric table too small: %dx%d", ...
                    fpath, size(M,1), size(M,2));
            end
            continue;
        end

        Cc = detect_cols_(M);
        id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
        if isempty(id)
            n_fail_match = n_fail_match + 1;
            if dbg_cfg.enable && numel(bad_match) < dbg_cfg.max_print
                bad_match(end+1) = sprintf("%s | no (iq,jq)=(%d,%d)", fpath, iq_pick, jq_pick);
            end
            continue;
        end

        Tvals(k) = round(H.T_K, 12);
        doping_header(k) = round(H.doping, 12);
        chiV(k) = M(id, Cc.Re);
        file(k) = fpath;

        n_ok = n_ok + 1;
    end

    mask = isfinite(Tvals) & isfinite(doping_header) & isfinite(chiV);

    S = struct();
    S.T      = Tvals(mask);
    S.doping = doping_header(mask);
    S.chi    = chiV(mask);
    S.file   = file(mask);

    dbg = struct();
    dbg.files_total = n;
    dbg.ok = n_ok;
    dbg.fail_header = n_fail_header;
    dbg.fail_read = n_fail_read;
    dbg.fail_cols = n_fail_cols;
    dbg.fail_match = n_fail_match;
    dbg.bad_header = bad_header;
    dbg.bad_read = bad_read;
    dbg.bad_cols = bad_cols;
    dbg.bad_match = bad_match;

    fprintf("[load scattered T-doping] files=%d ok=%d header_fail=%d read_fail=%d cols_fail=%d match_fail=%d\n", ...
        dbg.files_total, dbg.ok, dbg.fail_header, dbg.fail_read, dbg.fail_cols, dbg.fail_match);
end

function print_debug_summary_scattered_(dbg)
    fprintf("\n================ DEBUG SUMMARY (scattered T-doping) ================\n");
    fprintf("files_total=%d ok=%d header_fail=%d read_fail=%d cols_fail=%d match_fail=%d\n", ...
        dbg.files_total, dbg.ok, dbg.fail_header, dbg.fail_read, dbg.fail_cols, dbg.fail_match);
    fprintf("====================================================================\n\n");
end

% =====================================================================
% shared helpers
% =====================================================================
function H = parse_header_T_and_doping_robust_(fpath)
    H = struct('ok',false,'T_K',NaN,'doping',NaN);
    fid = fopen(fpath,'r');
    if fid < 0
        return;
    end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';
    for t = 1:4000
        ln = fgetl(fid);
        if ~ischar(ln)
            break;
        end
        s = strtrim(ln);
        if isempty(s)
            continue;
        end
        if s(1) ~= '#'
            break;
        end

        if ~isfinite(H.T_K)
            tok = regexp(s, ['^\#\s*T\s*=\s*' num], 'tokens','once');
            if isempty(tok)
                tok = regexp(s, ['^\#\s*T_K\s*=\s*' num], 'tokens','once');
            end
            if ~isempty(tok)
                H.T_K = str2double(tok{1});
            end
        end

        if ~isfinite(H.doping)
            tok = regexp(s, ['^\#\s*doping\s*=\s*' num], 'tokens','once');
            if ~isempty(tok)
                H.doping = str2double(tok{1});
            end
        end

        if isfinite(H.T_K) && isfinite(H.doping)
            break;
        end
    end

    H.ok = isfinite(H.T_K) && isfinite(H.doping);
end

function M = read_chi_numeric_skiphash_(fpath)
    fid = fopen(fpath,'r');
    if fid < 0
        error("Cannot open: %s", fpath);
    end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    rows = {};
    while true
        ln = fgetl(fid);
        if ~ischar(ln)
            break;
        end
        if isempty(ln)
            continue;
        end
        if ~isempty(regexp(ln, '^\s*#', 'once'))
            continue;
        end
        if isempty(strtrim(ln))
            continue;
        end
        v = sscanf(ln, '%f').';
        if isempty(v)
            continue;
        end
        rows{end+1,1} = v; %#ok<AGROW>
    end

    if isempty(rows)
        M = [];
        return;
    end

    ncol = max(cellfun(@numel, rows));
    M = nan(numel(rows), ncol);
    for i = 1:numel(rows)
        v = rows{i};
        M(i,1:numel(v)) = v;
    end

    lastFinite = find(any(isfinite(M),1), 1, 'last');
    if ~isempty(lastFinite)
        M = M(:,1:lastFinite);
    end
end

function C = detect_cols_(M)
    ncol = size(M,2);
    C = struct('iq',2,'jq',3,'Re',6);

    if ncol < 6
        C.iq = 1;
        C.jq = 2;
        C.Re = max(1, min(5, ncol));
        return;
    end

    intlike = @(x) all(abs(x - round(x)) < 1e-9);
    if ncol >= 7
        col1 = M(:,1);
        col2 = M(:,2);
        col3 = M(:,3);
        if intlike(col1) && (min(col1)==0 || min(col1)==1) && intlike(col2) && intlike(col3)
            if any(col2 < 0) || any(col3 < 0) || max(abs(col2)) > 5 || max(abs(col3)) > 5
                C.iq = 2;
                C.jq = 3;
                C.Re = 6;
                return;
            end
        end
    end

    C.iq = 1;
    C.jq = 2;
    C.Re = 5;
end