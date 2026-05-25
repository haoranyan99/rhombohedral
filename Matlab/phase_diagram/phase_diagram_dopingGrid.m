function out = phase_diagram_dopingGrid_fromMuHeader()
% phase_diagram_dopingGrid_fromMuHeader
%
% Goal:
%   Build phase diagram on a regular (doping, T) grid,
%   even if source files are organized by mu folders.
%
% Data source from each chi*.txt:
%   - T       : from header
%   - doping  : from header
%   - chi     : Re(chi) at chosen (iq,jq) from numeric table
%
% Then:
%   - evaluate coefficients using header doping
%   - classify phase 1/2/3/4
%   - interpolate / map onto a user-defined regular doping grid
%   - only plot within user-defined T and doping ranges
%
% Needs:
%   - make_realchi_params()
%   - make_realchi_coeff(par)

    clc; close all;

    % =========================================================
    % USER SETTINGS
    % =========================================================
    par  = make_realchi_params();
    coef = make_realchi_coeff(par);

    % q point to use
    iq_pick = round(par.iq_pick);
    jq_pick = round(par.jq_pick);

    % ---------- plot window ----------
    T_range   = [0, 10];         % [K]
    dop_range = [-2, 0];         % [10^12 cm^-2]

    % ---------- regular target grid ----------
    NT_grid   = 121;             % number of T grid points
    ND_grid   = 161;             % number of doping grid points

    % interpolation method for continuous fields
    interp_method = 'natural';   % 'nearest' | 'linear' | 'natural'
    % final phase assignment method on target grid:
    %   'nearest' : robust for categorical phase map
    %   'linear'  : interpolate phase index then round (less physical)
    phase_assign_method = 'nearest';

    % if multiple files land in same (T,doping) after rounding, average them
    T_round_digits   = 8;
    dop_round_digits = 6;

    % minimum number of valid points required
    min_valid_points = 5;

    % save figure
    save_png = true;

    % debug
    dbg = true;
    max_debug_print = 20;

    % =========================================================
    % SELECT ROOT
    % =========================================================
    default_root = "E:\rg_master\data\chi_400_1\D-0.012_multiT";
    if ~isfolder(default_root)
        default_root = string(pwd);
    end

    root = uigetdir(default_root, ...
        'Select root folder containing chi*.txt recursively');
    if isequal(root,0)
        error('User cancelled.');
    end
    root = string(root);

    fprintf('Root folder:\n  %s\n', root);
    fprintf('Use q-point: (iq,jq)=(%d,%d)\n', iq_pick, jq_pick);
    fprintf('Target T range      = [%.6g, %.6g] K\n', T_range(1), T_range(2));
    fprintf('Target doping range = [%.6g, %.6g]\n', dop_range(1), dop_range(2));

    % =========================================================
    % COLLECT SCATTERED POINTS
    % =========================================================
    files = dir(fullfile(root, "**", "chi*.txt"));
    if isempty(files)
        error('No chi*.txt found under:\n%s', root);
    end

    nfile = numel(files);

    T_raw    = nan(nfile,1);
    dop_raw  = nan(nfile,1);
    chi_raw  = nan(nfile,1);

    n_ok = 0;
    n_fail_header = 0;
    n_fail_read   = 0;
    n_fail_cols   = 0;
    n_fail_match  = 0;

    bad_header = strings(0);
    bad_read   = strings(0);
    bad_cols   = strings(0);
    bad_match  = strings(0);

    for k = 1:nfile
        fpath = fullfile(files(k).folder, files(k).name);

        % ---- parse header T + doping ----
        H = parse_header_T_and_doping_robust_(fpath);
        if ~H.ok
            n_fail_header = n_fail_header + 1;
            if dbg && numel(bad_header) < max_debug_print
                bad_header(end+1) = string(fpath);
            end
            continue;
        end

        % ---- read numeric table ----
        try
            M = read_chi_numeric_skiphash_(fpath);
        catch
            n_fail_read = n_fail_read + 1;
            if dbg && numel(bad_read) < max_debug_print
                bad_read(end+1) = string(fpath);
            end
            continue;
        end

        if isempty(M) || size(M,2) < 6
            n_fail_cols = n_fail_cols + 1;
            if dbg && numel(bad_cols) < max_debug_print
                bad_cols(end+1) = string(fpath);
            end
            continue;
        end

        % ---- detect columns ----
        try
            Cc = detect_cols_(M);
        catch
            n_fail_cols = n_fail_cols + 1;
            if dbg && numel(bad_cols) < max_debug_print
                bad_cols(end+1) = string(fpath);
            end
            continue;
        end

        % ---- find chosen q ----
        id = find(round(M(:,Cc.iq)) == iq_pick & round(M(:,Cc.jq)) == jq_pick, 1);
        if isempty(id)
            n_fail_match = n_fail_match + 1;
            if dbg && numel(bad_match) < max_debug_print
                bad_match(end+1) = string(fpath);
            end
            continue;
        end

        ReChi = M(id, Cc.Re);

        n_ok = n_ok + 1;
        T_raw(n_ok)   = round(H.T_K, T_round_digits);
        dop_raw(n_ok) = round(H.doping, dop_round_digits);
        chi_raw(n_ok) = ReChi;
    end

    T_raw   = T_raw(1:n_ok);
    dop_raw = dop_raw(1:n_ok);
    chi_raw = chi_raw(1:n_ok);

    fprintf('[collect] files=%d ok=%d | header_fail=%d read_fail=%d cols_fail=%d match_fail=%d\n', ...
        nfile, n_ok, n_fail_header, n_fail_read, n_fail_cols, n_fail_match);

    if dbg
        if ~isempty(bad_header)
            fprintf('\n[header_fail examples]\n');
            disp(bad_header);
        end
        if ~isempty(bad_match)
            fprintf('\n[match_fail examples]\n');
            disp(bad_match);
        end
    end

    if n_ok < min_valid_points
        error('Too few valid points (%d).', n_ok);
    end

    % =========================================================
    % KEEP ONLY USER-DEFINED WINDOW
    % =========================================================
    mask_win = ...
        T_raw   >= T_range(1)   & T_raw   <= T_range(2) & ...
        dop_raw >= dop_range(1) & dop_raw <= dop_range(2) & ...
        isfinite(chi_raw);

    T_raw   = T_raw(mask_win);
    dop_raw = dop_raw(mask_win);
    chi_raw = chi_raw(mask_win);

    if numel(T_raw) < min_valid_points
        error('Too few valid points inside requested T/doping window.');
    end

    fprintf('[window] kept=%d points in requested range.\n', numel(T_raw));

    % =========================================================
    % MERGE DUPLICATES IN (T,doping) BY AVERAGING chi
    % =========================================================
    TD = [T_raw, dop_raw];
    [TDu, ~, ic] = unique(TD, 'rows');

    T_sc   = TDu(:,1);
    dop_sc = TDu(:,2);
    chi_sc = accumarray(ic, chi_raw, [], @mean);

    % =========================================================
    % COMPUTE PHASE ON SCATTERED POINTS
    % =========================================================
    npt = numel(T_sc);

    phase_sc = zeros(npt,1);
    a1_sc    = nan(npt,1);
    a2_sc    = nan(npt,1);
    b2_sc    = nan(npt,1);
    crit_sc  = nan(npt,1);

    for i = 1:npt
        T0   = T_sc(i);
        dop0 = dop_sc(i);
        chi0 = chi_sc(i);

        try
            C = coef.eval(T0, dop0, chi0);
            a1 = C.a1;
            a2 = C.a2;
            b2 = C.b2;
        catch
            phase_sc(i) = 0;
            continue;
        end

        a2crit = (5/6) * b2^2;

        a1_sc(i)   = a1;
        a2_sc(i)   = a2;
        b2_sc(i)   = b2;
        crit_sc(i) = a2crit;

        if ~isfinite(a1) || ~isfinite(a2) || ~isfinite(b2) || ~isfinite(a2crit)
            phase_sc(i) = 0;
            continue;
        end

        if a2 >= a2crit
            if a1 < 0
                phase_sc(i) = 4;
            else
                phase_sc(i) = 1;
            end
        else
            if a1 < 0
                phase_sc(i) = 3;
            else
                phase_sc(i) = 2;
            end
        end
    end

    valid_sc = phase_sc > 0 & isfinite(T_sc) & isfinite(dop_sc);
    T_sc     = T_sc(valid_sc);
    dop_sc   = dop_sc(valid_sc);
    chi_sc   = chi_sc(valid_sc);
    phase_sc = phase_sc(valid_sc);
    a1_sc    = a1_sc(valid_sc);
    a2_sc    = a2_sc(valid_sc);
    b2_sc    = b2_sc(valid_sc);
    crit_sc  = crit_sc(valid_sc);

    if numel(T_sc) < min_valid_points
        error('Too few valid classified points after coefficient evaluation.');
    end

    % =========================================================
    % BUILD REGULAR TARGET GRID
    % =========================================================
    T_grid   = linspace(T_range(1),   T_range(2),   NT_grid);
    dop_grid = linspace(dop_range(1), dop_range(2), ND_grid);
    [DOP, TT] = meshgrid(dop_grid, T_grid);

    % continuous fields
    F_chi  = scatteredInterpolant(dop_sc, T_sc, chi_sc,  interp_method, 'none');
    F_a1   = scatteredInterpolant(dop_sc, T_sc, a1_sc,   interp_method, 'none');
    F_a2   = scatteredInterpolant(dop_sc, T_sc, a2_sc,   interp_method, 'none');
    F_b2   = scatteredInterpolant(dop_sc, T_sc, b2_sc,   interp_method, 'none');
    F_crit = scatteredInterpolant(dop_sc, T_sc, crit_sc, interp_method, 'none');

    chi_map  = F_chi(DOP, TT);
    a1_map   = F_a1(DOP, TT);
    a2_map   = F_a2(DOP, TT);
    b2_map   = F_b2(DOP, TT);
    crit_map = F_crit(DOP, TT);

    % categorical phase map
    switch lower(phase_assign_method)
        case 'nearest'
            F_phase = scatteredInterpolant(dop_sc, T_sc, double(phase_sc), 'nearest', 'none');
            phase_map = round(F_phase(DOP, TT));
        case 'linear'
            F_phase = scatteredInterpolant(dop_sc, T_sc, double(phase_sc), 'linear', 'none');
            phase_map = round(F_phase(DOP, TT));
        otherwise
            error('Unknown phase_assign_method: %s', phase_assign_method);
    end

    % remove points outside convex hull / undefined regions
    invalid = ~isfinite(chi_map) | ~isfinite(a1_map) | ~isfinite(a2_map) | ~isfinite(crit_map);
    phase_map(invalid) = 0;

    % =========================================================
    % PLOT
    % =========================================================
    fig = figure('Color','w', 'Name','Phase diagram on doping grid');
    ax = axes(fig); %#ok<LAXES>

    % colors: invalid + phase1..4
    c_invalid = [0.25 0.08 0.35];
    c1 = [0.20 0.55 0.95];
    c2 = [0.25 0.85 0.35];
    c3 = [1.00 0.60 0.10];
    c4 = [0.85 0.20 0.15];
    cmap = [c_invalid; c1; c2; c3; c4];

    dop_edges = centers_to_edges_(dop_grid(:));
    T_edges   = centers_to_edges_(T_grid(:));

    Cgrid = nan(numel(T_edges), numel(dop_edges));
    Cgrid(1:end-1,1:end-1) = phase_map;

    [DE, TE] = meshgrid(dop_edges, T_edges);
    surface(ax, DE, TE, zeros(size(DE)), Cgrid, ...
        'EdgeColor','none', 'FaceColor','flat');
    view(ax,2);
    set(ax,'YDir','normal');
    axis(ax,'tight');

    colormap(ax, cmap);
    caxis(ax, [0 4]);

    xlabel(ax, '$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$', 'Interpreter','latex');
    ylabel(ax, '$T\ (\mathrm{K})$', 'Interpreter','latex');
    title(ax, sprintf(['Phase diagram on regular doping grid | ' ...
        'q=(%d,%d) | source header doping, not folder mu'], iq_pick, jq_pick), ...
        'Interpreter','latex', 'FontWeight','normal');

    set(ax, 'FontSize', par.fontSize, ...
        'TickLabelInterpreter','latex', ...
        'LineWidth',1, 'TickDir','out', 'Box','on');

    % boundary contour
    hold(ax,'on');
    D = a2_map - crit_map;
    D(~isfinite(D)) = NaN;
    if any(isfinite(D(:)))
        contour(ax, dop_grid, T_grid, D, [0 0], 'k-', 'LineWidth', 1.2);
    end

    % overlay raw scattered points (optional, helpful for checking coverage)
    plot(ax, dop_sc, T_sc, 'k.', 'MarkerSize', 6);

    % legend
    h1 = plot(ax, nan,nan,'s','MarkerFaceColor',c1,'MarkerEdgeColor','none','MarkerSize',10);
    h2 = plot(ax, nan,nan,'s','MarkerFaceColor',c2,'MarkerEdgeColor','none','MarkerSize',10);
    h3 = plot(ax, nan,nan,'s','MarkerFaceColor',c3,'MarkerEdgeColor','none','MarkerSize',10);
    h4 = plot(ax, nan,nan,'s','MarkerFaceColor',c4,'MarkerEdgeColor','none','MarkerSize',10);
    h0 = plot(ax, nan,nan,'s','MarkerFaceColor',c_invalid,'MarkerEdgeColor','none','MarkerSize',10);
    hs = plot(ax, nan,nan,'k.', 'MarkerSize', 10);
    hold(ax,'off');

    legend([h1 h2 h3 h4 h0 hs], ...
        {'phase1', 'phase2', 'phase3', 'phase4', 'invalid', 'raw points'}, ...
        'Interpreter','latex', 'Location','eastoutside');

    % datatip
    dcm = datacursormode(fig);
    set(dcm, 'Enable','on', 'SnapToDataVertex','off', 'DisplayStyle','datatip');
    set(dcm, 'UpdateFcn', @tip_cb_);

    % =========================================================
    % SAVE
    % =========================================================
    out_png = "";
    if save_png
        out_dir = fullfile(root, "..", "plot");
        if ~exist(out_dir, "dir"), mkdir(out_dir); end
        out_png = fullfile(out_dir, ...
            sprintf('phase_dopingGrid_iq%d_jq%d_T%.3g_%.3g_dop%.3g_%.3g.png', ...
            iq_pick, jq_pick, T_range(1), T_range(2), dop_range(1), dop_range(2)));
        exportgraphics(fig, out_png, 'Resolution', 300);
        fprintf('Saved: %s\n', out_png);
    end

    % =========================================================
    % OUTPUT
    % =========================================================
    out = struct();
    out.root      = root;
    out.iq_pick   = iq_pick;
    out.jq_pick   = jq_pick;
    out.T_grid    = T_grid;
    out.dop_grid  = dop_grid;
    out.phase_map = phase_map;
    out.a1_map    = a1_map;
    out.a2_map    = a2_map;
    out.b2_map    = b2_map;
    out.crit_map  = crit_map;
    out.chi_map   = chi_map;
    out.T_sc      = T_sc;
    out.dop_sc    = dop_sc;
    out.phase_sc  = phase_sc;
    out.png       = out_png;

    % =========================================================
    % DATATIP CALLBACK
    % =========================================================
    function txt = tip_cb_(~, evt)
        pos = evt.Position;
        x = pos(1); % doping
        y = pos(2); % T

        iu = find(x >= dop_edges(1:end-1) & x < dop_edges(2:end), 1, 'first');
        it = find(y >= T_edges(1:end-1) & y < T_edges(2:end), 1, 'first');

        if isempty(iu) || isempty(it)
            txt = {'(out of range)'};
            return;
        end

        ph  = phase_map(it, iu);
        T0  = T_grid(it);
        d0  = dop_grid(iu);
        a10 = a1_map(it, iu);
        a20 = a2_map(it, iu);
        b20 = b2_map(it, iu);
        cr0 = crit_map(it, iu);
        ch0 = chi_map(it, iu);

        txt = {
            sprintf('T = %.6g K', T0)
            sprintf('doping = %.6g', d0)
            sprintf('phase = %d (%s)', ph, phase_name1234_(ph))
            sprintf('Re(chi) = %.6g', ch0)
            sprintf('a1 = %.6g', a10)
            sprintf('a2 = %.6g', a20)
            sprintf('b2 = %.6g', b20)
            sprintf('a2crit = %.6g', cr0)
        };
    end
end

% =====================================================================
% helpers
% =====================================================================
function edges = centers_to_edges_(c)
    c = c(:);
    if numel(c) == 1
        dc = 1;
        edges = [c(1)-dc/2; c(1)+dc/2];
        return;
    end
    mid = (c(1:end-1) + c(2:end))/2;
    edges = [c(1) - (mid(1)-c(1)); mid; c(end) + (c(end)-mid(end))];
end

function name = phase_name1234_(ph)
    switch ph
        case 1, name = 'phase1';
        case 2, name = 'phase2';
        case 3, name = 'phase3';
        case 4, name = 'phase4';
        otherwise, name = 'invalid/missing';
    end
end

function H = parse_header_T_and_doping_robust_(fpath)
    H = struct('ok',false,'T_K',NaN,'doping',NaN);

    fid = fopen(fpath,'r');
    if fid < 0, return; end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';

    for t = 1:4000
        ln = fgetl(fid);
        if ~ischar(ln), break; end

        s = strtrim(ln);
        if isempty(s), continue; end
        if s(1) ~= '#', break; end

        if ~isfinite(H.T_K)
            tok = regexp(s, ['^\#\s*T\s*=\s*' num '\b'], 'tokens', 'once');
            if isempty(tok)
                tok = regexp(s, ['^\#\s*T_K\s*=\s*' num '\b'], 'tokens', 'once');
            end
            if isempty(tok)
                tok = regexp(s, ['^\#\s*T\s*=\s*' num], 'tokens', 'once');
            end
            if ~isempty(tok)
                H.T_K = str2double(tok{1});
            end
        end

        if ~isfinite(H.doping)
            tok = regexp(s, ['^\#\s*doping\s*=\s*' num '\b'], 'tokens', 'once', 'ignorecase');
            if isempty(tok)
                tok = regexp(s, ['^\#\s*dop\s*=\s*' num '\b'], 'tokens', 'once', 'ignorecase');
            end
            if isempty(tok)
                tok = regexp(s, ['^\#\s*doping\s*=\s*' num], 'tokens', 'once', 'ignorecase');
            end
            if isempty(tok)
                tok = regexp(s, ['^\#\s*dop\s*=\s*' num], 'tokens', 'once', 'ignorecase');
            end
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
    if fid < 0, error("Cannot open: %s", fpath); end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    rows = {};
    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        if isempty(ln), continue; end
        if ~isempty(regexp(ln, '^\s*#', 'once')), continue; end
        if isempty(strtrim(ln)), continue; end
        v = sscanf(ln, '%f').';
        if isempty(v), continue; end
        rows{end+1,1} = v; %#ok<AGROW>
    end

    if isempty(rows), M = []; return; end
    ncol = max(cellfun(@numel, rows));
    M = nan(numel(rows), ncol);
    for i = 1:numel(rows)
        v = rows{i};
        M(i,1:numel(v)) = v;
    end
    lastFinite = find(any(isfinite(M),1), 1, 'last');
    if ~isempty(lastFinite), M = M(:,1:lastFinite); end
end

function C = detect_cols_(M)
    ncol = size(M,2);

    if ncol < 6
        error('detect_cols_: numeric table has fewer than 6 columns.');
    end

    intlike = @(x) all(abs(x - round(x)) < 1e-9);

    col1 = M(:,1);
    col2 = M(:,2);
    col3 = M(:,3);

    % Common format:
    % idx iq jq qx qy Re Im ...
    if intlike(col1) && intlike(col2) && intlike(col3)
        C.iq = 2;
        C.jq = 3;
        C.Re = 6;
        return;
    end

    % Compact format:
    % iq jq qx qy Re Im
    if intlike(col1) && intlike(col2)
        C.iq = 1;
        C.jq = 2;
        C.Re = 5;
        return;
    end

    % fallback
    C.iq = 2;
    C.jq = 3;
    C.Re = 6;
end