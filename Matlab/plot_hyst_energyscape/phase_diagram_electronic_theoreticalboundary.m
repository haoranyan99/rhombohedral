function out = phase_diagram_electronic_theoreticalboundary()
% phase_diagram_electronic_doping_nn
% Electronic phase diagram using doping as x-axis.
%
% Phase definition:
%   a1 > 0  -> red phase
%   a1 < 0  -> blue phase
%
% Data source:
%   read chi*.txt recursively
%
% Important:
%   original data are aligned in mu-grid for each T, but not in doping-grid.
%   therefore this script builds a common doping axis and renders by nearest
%   neighbor assignment row-by-row.

    % =========================================================
    % 1. params / coeff / style
    % =========================================================
    par  = make_realchi_params();
    coef = make_realchi_coeff(par);

    fs = 18;

    % user controls for rendered doping grid
    dop_min_user = -4.0;
    dop_max_user =  0.0;
    Ndop_render  = 401;

    % nearest-neighbor tolerance factor:
    % if nearest raw point is farther than tol_factor * local grid spacing,
    % leave that rendered point as NaN.
    tol_factor = 0.75;

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
    % 3. load scattered data from files
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
    a1_raw    = nan(Nraw,1);
    phase_raw = nan(Nraw,1);   % 1=red (a1>0), 2=blue (a1<0)

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

        if a1 > 0
            phase_raw(k) = 1;   % red
        elseif a1 < 0
            phase_raw(k) = 2;   % blue
        else
            phase_raw(k) = NaN; % exactly zero, leave blank
        end
    end

    valid_mask = isfinite(S.T) & isfinite(S.doping) & isfinite(a1_raw) & isfinite(phase_raw);
    T_raw      = S.T(valid_mask);
    dop_raw    = S.doping(valid_mask);
    chi_raw    = S.chi(valid_mask);
    a1_raw     = a1_raw(valid_mask);
    phase_raw  = phase_raw(valid_mask);
    file_raw   = S.file(valid_mask);

    if isempty(T_raw)
        error('All raw points became invalid after evaluating a1.');
    end

    % =========================================================
    % 5. build rendered regular grid in (T, doping)
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

    phase_map = nan(NT, Ndop);
    a1_map    = nan(NT, Ndop);
    chi_map   = nan(NT, Ndop);

    % also keep nearest raw doping and file for datatip
    dop_src_map  = nan(NT, Ndop);
    dist_map     = nan(NT, Ndop);
    file_map     = strings(NT, Ndop);

    for iT = 1:NT
        T0 = T_list(iT);
        mT = (T_raw == T0);

        dop_row   = dop_raw(mT);
        a1_row    = a1_raw(mT);
        chi_row   = chi_raw(mT);
        phase_row = phase_raw(mT);
        file_row  = file_raw(mT);

        if isempty(dop_row)
            continue;
        end

        % sort by doping
        [dop_row, order] = sort(dop_row);
        a1_row    = a1_row(order);
        chi_row   = chi_row(order);
        phase_row = phase_row(order);
        file_row  = file_row(order);

        % remove duplicate doping points by averaging a1/chi and majority-like phase
        [dop_u, ~, ic] = unique(round(dop_row, 10));
        a1_u    = nan(size(dop_u));
        chi_u   = nan(size(dop_u));
        phase_u = nan(size(dop_u));
        file_u  = strings(size(dop_u));

        for j = 1:numel(dop_u)
            mm = (ic == j);
            a1_u(j)  = mean(a1_row(mm), 'omitnan');
            chi_u(j) = mean(chi_row(mm), 'omitnan');

            phs = phase_row(mm);
            phs = phs(isfinite(phs));
            if isempty(phs)
                phase_u(j) = NaN;
            else
                phase_u(j) = round(mean(phs)); % should still be 1 or 2
            end

            idx_first = find(mm, 1, 'first');
            if ~isempty(idx_first)
                file_u(j) = file_row(idx_first);
            end
        end

        dop_u = double(dop_u);

        % local tolerance based on rendered spacing and raw spacing
        if numel(dop_u) >= 2
            raw_spacing = median(diff(dop_u));
        else
            raw_spacing = inf;
        end
        render_spacing = mean(diff(dop_list));
        tol = tol_factor * max(raw_spacing, render_spacing);

        % nearest-neighbor render
        for j = 1:Ndop
            x = dop_list(j);
            [dmin, idx] = min(abs(dop_u - x));

            if isempty(idx) || ~isfinite(dmin)
                continue;
            end

            if dmin <= tol
                phase_map(iT, j)   = phase_u(idx);
                a1_map(iT, j)      = a1_u(idx);
                chi_map(iT, j)     = chi_u(idx);
                dop_src_map(iT, j) = dop_u(idx);
                dist_map(iT, j)    = dmin;
                file_map(iT, j)    = file_u(idx);
            end
        end
    end

    fprintf('Rendered valid points = %d / %d\n', nnz(isfinite(phase_map)), numel(phase_map));

    % =========================================================
    % 6. plot
    % =========================================================
    fig = figure('Color','w', 'Name','Electronic phase diagram (doping axis)');
    ax = axes(fig);
    hold(ax, 'on');

    img_handle = [];

    Cplot = nan(size(phase_map));
    Cplot(phase_map == 1) = 1;   % red
    Cplot(phase_map == 2) = 2;   % blue

    img_handle = imagesc(ax, dop_list, T_list, Cplot);
    set(ax, 'YDir', 'normal');

    colormap(ax, [0.85 0.20 0.15;
                  0.20 0.45 0.90]);
    caxis(ax, [1 2]);

    cb = colorbar(ax);
    cb.Ticks = [1.25, 1.75];
    cb.TickLabels = {'off', 'on'};
    cb.FontSize = fs;

    xlabel(ax, 'doping (10^{12} cm^{-2})', 'FontSize', fs);
    ylabel(ax, 'T (K)', 'FontSize', fs);
    set(ax, 'FontSize', fs, 'LineWidth', 1.2, 'Box', 'on');
    axis(ax, 'tight');

    title(ax, 'Electronic phase transition', 'FontSize', fs, 'FontWeight', 'normal');

    % =========================================================
    % 7. datatip
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
    % 8. output
    % =========================================================
    out = struct();
    out.dop_list    = dop_list;
    out.T_list      = T_list;
    out.phase_map   = phase_map;
    out.a1_map      = a1_map;
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

        dop0 = dop_list(j);
        T0   = T_list(i);

        a10   = a1_map(i,j);
        chi0  = chi_map(i,j);
        ph    = phase_map(i,j);
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

        if ph == 1
            phname = 'red (a1 > 0)';
        else
            phname = 'blue (a1 < 0)';
        end

        txt = {
            sprintf('doping(render) = %.6g', dop0)
            sprintf('T = %.6g', T0)
            sprintf('phase = %s', phname)
            sprintf('a1 = %.6g', a10)
            sprintf('chi = %.6g', chi0)
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