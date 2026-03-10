function out = phase_diagram_LatticePhaseTransition()
% phase_diagram_group14_vs_23_realchi
% Grouped phase map:
%   groupA = phases {1,4}  (same color)
%   groupB = phases {2,3}  (same color)
%   invalid/missing = 0
%
% Classification (same as your 1-4 logic):
%   a2crit = (5/6) * b2^2
%   if a2 >= a2crit: phase-family {1 or 4} by sign(a1)
%   else            phase-family {2 or 3} by sign(a1)
%
% Plot:
%   - aligned blocks by true coordinate edges (surface w/ flat face)
%   - ONLY contour for boundary a2-a2crit=0 (group 14 vs 23 boundary)
%   - axis tick modes = auto
%   - datatip shows (T, folder U, group, phase, header doping, a1,a2,b2,a2crit)
%
% Needs:
%   - make_realchi_params()
%   - make_realchi_coeff(par): coef.eval(T, doping, chi)->a1,b1,a2,b2,c2
%   - load_chi_grid_folderUI_debug_() and its helper subfunctions
%     (copied from your previous script)

    % -------------------- params / coeff --------------------
    par  = make_realchi_params();
    coef = make_realchi_coeff(par);

    % ---------------- choose folder -----------------
    default_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_200\D0.067";
    if ~isfolder(default_root)
        warning("Default folder not found: %s\nFallback to pwd.", default_root);
        default_root = string(pwd);
    end

    root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
    if isequal(root, 0), error('User cancelled.'); end
    root = string(root);
    fprintf("Root folder:\n  %s\n", root);

    iq_pick = round(par.iq_pick);
    jq_pick = round(par.jq_pick);
    fprintf("Pick absolute q-point: (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);

    % ---------------- debug config ----------------
    dbg = struct();
    dbg.enable = true;
    dbg.max_print_per_reason = 80;
    dbg.max_missing_print = 40;

    % ---------------- load grid (folder UI logic) ----------------
    [G, dbg_out] = load_chi_grid_folderUI_debug_(root, iq_pick, jq_pick, dbg);

    T_list  = G.T_list;           % header T
    U_list  = G.U_list;           % folder mu/doping
    u_tag   = G.u_tag;            % 'mu' or 'doping'
    chi_map = G.chi_map;          % Re(chi)
    dop_map = G.doping_map;       % header doping (for coef.eval)

    NT = numel(T_list);
    NU = numel(U_list);

    fprintf("Loaded grid: NT=%d, N%s=%d, valid=%d | q=(%d,%d)\n", ...
        NT, u_tag, NU, nnz(isfinite(chi_map) & isfinite(dop_map)), iq_pick, jq_pick);

    if dbg.enable
        print_debug_summary_folderUI_(dbg_out, dbg, G);
    end

    if NT==0 || NU==0
        error("Parsed 0 valid points. Check folder tree / headers / (iq,jq).");
    end

    % ---------------- classify 1..4 then group into 0/1/2 ----------------
    % phase code:
    %   0 invalid
    %   1 phase1
    %   2 phase2
    %   3 phase3
    %   4 phase4
    phase_map = zeros(NT, NU);

    a1_map   = nan(NT, NU);
    a2_map   = nan(NT, NU);
    b2_map   = nan(NT, NU);
    crit_map = nan(NT, NU);

    % group_map:
    %   0 invalid
    %   1 group {1,4}
    %   2 group {2,3}
    group_map = zeros(NT, NU);

    % optional parallel
    try
        p = gcp('nocreate');
        if isempty(p), parpool; end
        usePar = true;
    catch
        fprintf("[warn] Parallel pool not available. Running serial.\n");
        usePar = false;
    end

    if usePar
        parfor iT = 1:NT
            T = T_list(iT);

            ph_row  = zeros(1,NU);
            g_row   = zeros(1,NU);

            a1_row  = nan(1,NU);
            a2_row  = nan(1,NU);
            b2_row  = nan(1,NU);
            cr_row  = nan(1,NU);

            for iU = 1:NU
                chi_used = chi_map(iT,iU);
                dop = dop_map(iT,iU);

                if ~isfinite(chi_used) || ~isfinite(dop)
                    ph_row(iU) = 0;
                    g_row(iU)  = 0;
                    continue;
                end

                try
                    C  = coef.eval(T, dop, chi_used);
                    a1 = C.a1; a2 = C.a2; b2 = C.b2;
                catch
                    ph_row(iU) = 0;
                    g_row(iU)  = 0;
                    continue;
                end

                a2crit = (5/6) * (b2^2);

                a1_row(iU) = a1;
                a2_row(iU) = a2;
                b2_row(iU) = b2;
                cr_row(iU) = a2crit;

                if ~isfinite(a1) || ~isfinite(a2) || ~isfinite(b2) || ~isfinite(a2crit)
                    ph = 0;
                else
                    if a2 >= a2crit
                        ph = 1;      % phase1
                        if a1 < 0, ph = 4; end
                    else
                        ph = 2;      % phase2
                        if a1 < 0, ph = 3; end
                    end
                end

                ph_row(iU) = ph;

                % group
                if ph==1 || ph==4
                    g_row(iU) = 1;
                elseif ph==2 || ph==3
                    g_row(iU) = 2;
                else
                    g_row(iU) = 0;
                end
            end

            phase_map(iT,:) = ph_row;
            group_map(iT,:) = g_row;

            a1_map(iT,:)    = a1_row;
            a2_map(iT,:)    = a2_row;
            b2_map(iT,:)    = b2_row;
            crit_map(iT,:)  = cr_row;
        end
    else
        for iT = 1:NT
            T = T_list(iT);
            for iU = 1:NU
                chi_used = chi_map(iT,iU);
                dop = dop_map(iT,iU);

                if ~isfinite(chi_used) || ~isfinite(dop)
                    phase_map(iT,iU) = 0;
                    group_map(iT,iU) = 0;
                    continue;
                end

                try
                    C  = coef.eval(T, dop, chi_used);
                    a1 = C.a1; a2 = C.a2; b2 = C.b2;
                catch
                    phase_map(iT,iU) = 0;
                    group_map(iT,iU) = 0;
                    continue;
                end

                a2crit = (5/6) * (b2^2);

                a1_map(iT,iU)   = a1;
                a2_map(iT,iU)   = a2;
                b2_map(iT,iU)   = b2;
                crit_map(iT,iU) = a2crit;

                if ~isfinite(a1) || ~isfinite(a2) || ~isfinite(b2) || ~isfinite(a2crit)
                    ph = 0;
                else
                    if a2 >= a2crit
                        ph = 1;
                        if a1 < 0, ph = 4; end
                    else
                        ph = 2;
                        if a1 < 0, ph = 3; end
                    end
                end
                phase_map(iT,iU) = ph;

                if ph==1 || ph==4
                    group_map(iT,iU) = 1;
                elseif ph==2 || ph==3
                    group_map(iT,iU) = 2;
                else
                    group_map(iT,iU) = 0;
                end
            end
        end
    end

    % ---------------- plot: group map + ONLY boundary (a2-a2crit=0) ----------------
    fig = figure('Color','w','Name','Phase groups: {1,4} vs {2,3}');
    ax = axes(fig); %#ok<LAXES>

    % --- colors ---
    c_invalid = [0.25 0.08 0.35];

    % 你要的：phase14 红、phase23 蓝（可随时对调）
    c14 = [0.85 0.20 0.15];   % group {1,4}
    c23 = [0.20 0.55 0.95];   % group {2,3}

    cmap2 = [c_invalid; c14; c23];  % values 0,1,2

    U_edges = centers_to_edges_(U_list);
    T_edges = centers_to_edges_(T_list);

    Cgrid = nan(numel(T_edges), numel(U_edges));
    Cgrid(1:end-1, 1:end-1) = group_map;

    [UUe, TTe] = meshgrid(U_edges, T_edges);
    surface(ax, UUe, TTe, zeros(size(UUe)), Cgrid, 'EdgeColor','none', 'FaceColor','flat');
    view(ax, 2);
    axis(ax, 'tight');
    set(ax,'YDir','normal');

    colormap(ax, cmap2);
    caxis(ax, [0 2]);

    ax.XTickMode = 'auto';
    ax.YTickMode = 'auto';

    if strcmpi(u_tag,'doping')
        xlabel(ax, '$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$','Interpreter','latex');
    else
        xlabel(ax, '$\mu\ (\mathrm{eV})$','Interpreter','latex');
    end
    ylabel(ax, '$T$','Interpreter','latex');

    title(ax, sprintf('Grouped phases: {1,4} vs {2,3} | boundary: $a_2-\\frac{5}{6}b_2^2=0$ | q=(%d,%d) | x=%s(folder)', ...
        iq_pick, jq_pick, u_tag), 'Interpreter','latex','FontWeight','normal');

    set(ax,'FontSize',par.fontSize,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out','Box','on');

    % --- boundary contour (a2-a2crit=0) ---
    hold(ax,'on');
    [UUc, TTc] = meshgrid(U_list, T_list);
    D = a2_map - crit_map;
    D(~isfinite(D)) = NaN;
    if any(isfinite(D(:)))
        contour(ax, UUc, TTc, D, [0 0], 'k-', 'LineWidth', 1.4);
    end
    hold(ax,'off');

    % --- legend (no colorbar) ---
    hold(ax,'on');
    h14 = plot(ax, nan,nan,'s','MarkerFaceColor',c14,'MarkerEdgeColor','none','MarkerSize',10);
    h23 = plot(ax, nan,nan,'s','MarkerFaceColor',c23,'MarkerEdgeColor','none','MarkerSize',10);
    h0  = plot(ax, nan,nan,'s','MarkerFaceColor',c_invalid,'MarkerEdgeColor','none','MarkerSize',10);
    hold(ax,'off');

    legend([h14 h23 h0], ...
        {'group A: phase1 or phase4 ($a_2 \ge \frac{5}{6}b_2^2$)', ...
         'group B: phase2 or phase3 ($a_2 <  \frac{5}{6}b_2^2$)', ...
         'invalid / missing'}, ...
         'Interpreter','latex','Location','eastoutside');

    % --- datatip / hover ---
    dcm = datacursormode(fig);
    set(dcm, 'Enable','on', 'SnapToDataVertex','off', 'DisplayStyle','datatip');
    set(dcm, 'UpdateFcn', @tip_cb_);

    % ---------------- save ----------------
    out_dir = fullfile(root, "..", "plot");
    if ~exist(out_dir,"dir"), mkdir(out_dir); end
    out_png = fullfile(out_dir, sprintf("phase_group14_vs_23_iq%d_jq%d_x%s.png", ...
        iq_pick, jq_pick, u_tag));
    exportgraphics(fig, out_png, 'Resolution', 300);
    fprintf("Saved: %s\n", out_png);

    % output
    out = struct();
    out.group_map  = group_map;
    out.phase_map  = phase_map;
    out.T_list     = T_list;
    out.U_list     = U_list;
    out.u_tag      = u_tag;
    out.chi_map    = chi_map;
    out.doping_map = dop_map;
    out.a1_map     = a1_map;
    out.a2_map     = a2_map;
    out.b2_map     = b2_map;
    out.crit_map   = crit_map;
    out.iq_pick    = iq_pick;
    out.jq_pick    = jq_pick;
    out.png        = out_png;

    % ===================== nested: datatip callback =====================
    function txt = tip_cb_(~, evt)
        pos = evt.Position;
        x = pos(1); y = pos(2);

        iu = find(x >= U_edges(1:end-1) & x < U_edges(2:end), 1, 'first');
        it = find(y >= T_edges(1:end-1) & y < T_edges(2:end), 1, 'first');

        if isempty(iu) || isempty(it)
            txt = {'(out of range)'};
            return;
        end

        g    = group_map(it, iu);
        ph   = phase_map(it, iu);
        T0   = T_list(it);
        U0   = U_list(iu);
        dop0 = dop_map(it, iu);

        a10  = a1_map(it, iu);
        a20  = a2_map(it, iu);
        b20  = b2_map(it, iu);
        cr0  = crit_map(it, iu);

        if strcmpi(u_tag,'doping')
            Ustr = sprintf('doping(folder) = %.4f', U0);
        else
            Ustr = sprintf('mu(folder) = %.6g', U0);
        end

        txt = {
            sprintf('T = %.6g K', T0)
            Ustr
            sprintf('group = %d (%s)', g, group_name_(g))
            sprintf('phase = %d', ph)
            sprintf('doping(header) = %.6g', dop0)
            sprintf('a1 = %.6g', a10)
            sprintf('a2 = %.6g', a20)
            sprintf('b2 = %.6g', b20)
            sprintf('a2crit = (5/6)*b2^2 = %.6g', cr0)
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

function s = group_name_(g)
    switch g
        case 1, s = '{phase1, phase4}';
        case 2, s = '{phase2, phase3}';
        otherwise, s = 'invalid/missing';
    end
end

% =====================================================================
% Folder-UI loader (mu/doping from path), doping from header, chi from table
% (copied from your current working version)
% =====================================================================
function [G, dbg] = load_chi_grid_folderUI_debug_(root_dir, iq_pick, jq_pick, dbg)
    root_dir = string(root_dir);
    L = dir(fullfile(root_dir, "**", "chi*.txt"));
    if isempty(L), error("No chi*.txt found under: %s", root_dir); end

    iq_pick = round(iq_pick);
    jq_pick = round(jq_pick);

    n = numel(L);

    Tvals = nan(n,1);
    Uvals = nan(n,1);
    chiV  = nan(n,1);
    doping_header = nan(n,1);

    Ttag = strings(n,1);
    Utag = strings(n,1);

    u_tag = "";

    n_ok=0; n_fail_header=0; n_fail_read=0; n_fail_match=0; n_fail_cols=0; n_fail_path=0;
    bad_header=strings(0); bad_read=strings(0); bad_match=strings(0); bad_cols=strings(0); bad_path=strings(0);

    for k = 1:n
        fpath = fullfile(L(k).folder, L(k).name);

        [tTag, uTag, uKind, uVal] = parse_T_and_U_from_path_(string(fpath));
        Ttag(k) = tTag;
        Utag(k) = uTag;

        if u_tag == "" && uKind ~= ""
            u_tag = uKind;
        end

        if uKind == "" || ~isfinite(uVal)
            n_fail_path = n_fail_path + 1;
            if dbg.enable && numel(bad_path) < dbg.max_print_per_reason
                bad_path(end+1) = sprintf("%s | cannot parse folder U | %s/%s", string(fpath), tTag, uTag);
            end
            continue;
        end

        H = parse_header_T_and_doping_robust_(fpath);
        if ~H.ok
            n_fail_header = n_fail_header + 1;
            if dbg.enable && numel(bad_header) < dbg.max_print_per_reason
                bad_header(end+1) = sprintf("%s | header missing T/doping (T=%g doping=%g) | %s/%s", ...
                    string(fpath), H.T_K, H.doping, tTag, uTag);
            end
            continue;
        end

        try
            M = read_chi_numeric_skiphash_(fpath);
        catch ME
            n_fail_read = n_fail_read + 1;
            if dbg.enable && numel(bad_read) < dbg.max_print_per_reason
                bad_read(end+1) = sprintf("%s | read failed: %s | %s/%s", string(fpath), string(ME.message), tTag, uTag);
            end
            continue;
        end

        if isempty(M) || size(M,2) < 6
            n_fail_cols = n_fail_cols + 1;
            if dbg.enable && numel(bad_cols) < dbg.max_print_per_reason
                bad_cols(end+1) = sprintf("%s | numeric table too small: %dx%d | %s/%s", string(fpath), size(M,1), size(M,2), tTag, uTag);
            end
            continue;
        end

        Cc = detect_cols_(M);
        id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
        if isempty(id)
            n_fail_match = n_fail_match + 1;
            if dbg.enable && numel(bad_match) < dbg.max_print_per_reason
                bad_match(end+1) = sprintf("%s | no (iq,jq)=(%d,%d) | %s/%s", string(fpath), iq_pick, jq_pick, tTag, uTag);
            end
            continue;
        end

        ReChi = M(id, Cc.Re);

        Tvals(k) = round(H.T_K, 12);
        Uvals(k) = quantize_U_(uKind, uVal);
        chiV(k)  = ReChi;
        doping_header(k) = round(H.doping, 12);

        n_ok = n_ok + 1;
    end

    if u_tag == ""
        u_tag = "doping";
    end

    dbg.files_total = n;
    dbg.ok = n_ok;
    dbg.fail_path = n_fail_path;
    dbg.fail_header = n_fail_header;
    dbg.fail_read = n_fail_read;
    dbg.fail_cols = n_fail_cols;
    dbg.fail_match = n_fail_match;
    dbg.bad_path = bad_path;
    dbg.bad_header = bad_header;
    dbg.bad_read = bad_read;
    dbg.bad_cols = bad_cols;
    dbg.bad_match = bad_match;

    fprintf("[load folderUI] files=%d ok=%d path_fail=%d header_fail=%d read_fail=%d cols_fail=%d match_fail=%d | axis=%s\n", ...
        dbg.files_total, dbg.ok, dbg.fail_path, dbg.fail_header, dbg.fail_read, dbg.fail_cols, dbg.fail_match, u_tag);

    mask = isfinite(Tvals) & isfinite(Uvals) & isfinite(chiV) & isfinite(doping_header);
    Tvals = Tvals(mask);
    Uvals = Uvals(mask);
    chiV  = chiV(mask);
    doping_header = doping_header(mask);
    Ttag = Ttag(mask);
    Utag = Utag(mask);

    if isempty(Tvals)
        G = struct('T_list',zeros(0,1),'U_list',zeros(0,1), ...
                   'chi_map',[],'doping_map',[], ...
                   'u_tag',u_tag,'Tfolder_map',strings(0), 'Ufolder_map',strings(0));
        return;
    end

    T_list = sort(unique(Tvals));
    U_list = sort(unique(Uvals));
    NT = numel(T_list);
    NU = numel(U_list);

    chi_sum = zeros(NT,NU);
    chi_cnt = zeros(NT,NU);
    dop_sum = zeros(NT,NU);
    dop_cnt = zeros(NT,NU);

    Tfolder_map = strings(NT,NU);
    Ufolder_map = strings(NT,NU);

    for i = 1:numel(Tvals)
        iT = find(T_list==Tvals(i),1);
        iU = find(U_list==Uvals(i),1);
        if isempty(iT) || isempty(iU), continue; end

        chi_sum(iT,iU) = chi_sum(iT,iU) + chiV(i);
        chi_cnt(iT,iU) = chi_cnt(iT,iU) + 1;

        dop_sum(iT,iU) = dop_sum(iT,iU) + doping_header(i);
        dop_cnt(iT,iU) = dop_cnt(iT,iU) + 1;

        if Tfolder_map(iT,iU) == ""
            Tfolder_map(iT,iU) = Ttag(i);
            Ufolder_map(iT,iU) = Utag(i);
        end
    end

    chi_map = nan(NT,NU);
    m = chi_cnt>0;
    chi_map(m) = chi_sum(m)./chi_cnt(m);

    doping_map = nan(NT,NU);
    md = dop_cnt>0;
    doping_map(md) = dop_sum(md)./dop_cnt(md);

    G = struct('T_list',T_list,'U_list',U_list, ...
               'chi_map',chi_map,'doping_map',doping_map, ...
               'u_tag',u_tag,'Tfolder_map',Tfolder_map,'Ufolder_map',Ufolder_map);
end

function [Tfolder, Ufolder, ukind, uval] = parse_T_and_U_from_path_(fpath)
    Tfolder = ""; Ufolder = ""; ukind = ""; uval = NaN;
    p = replace(string(fpath), "\", "/");

    tokT = regexp(p, '/(T[^/]+)/', 'tokens', 'once');
    if ~isempty(tokT), Tfolder = string(tokT{1}); end

    extract_num = @(s) str2double(regexp(s, '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)', 'match', 'once'));

    tokMu = regexp(p, '/(mu[^/]+)/', 'tokens', 'once');
    if ~isempty(tokMu)
        Ufolder = string(tokMu{1});
        ukind = "mu";
        uval = extract_num(Ufolder);
        if isfinite(uval), return; end
        Ufolder=""; ukind=""; uval=NaN;
    end

    tokDp = regexp(p, '/(doping[^/]+)/', 'tokens', 'once');
    if ~isempty(tokDp)
        Ufolder = string(tokDp{1});
        ukind = "doping";
        uval = extract_num(Ufolder);
        if isfinite(uval), return; end
        Ufolder=""; ukind=""; uval=NaN;
    end
end

function u = quantize_U_(ukind, uval)
    if ukind == "doping"
        u = round(uval, 4);
    else
        u = round(uval, 12);
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
            tok = regexp(s, ['^\#\s*T\s*=\s*' num], 'tokens','once');
            if isempty(tok)
                tok = regexp(s, ['^\#\s*T_K\s*=\s*' num], 'tokens','once');
            end
            if ~isempty(tok), H.T_K = str2double(tok{1}); end
        end

        if ~isfinite(H.doping)
            tok = regexp(s, ['^\#\s*doping\s*=\s*' num], 'tokens','once');
            if ~isempty(tok), H.doping = str2double(tok{1}); end
        end

        if isfinite(H.T_K) && isfinite(H.doping), break; end
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
    C = struct('iq',2,'jq',3,'Re',6);
    if ncol < 6
        C.iq = 1; C.jq = 2; C.Re = max(1, min(5,ncol));
        return;
    end
    intlike = @(x) all(abs(x - round(x)) < 1e-9);
    if ncol >= 7
        col1 = M(:,1); col2 = M(:,2); col3 = M(:,3);
        if intlike(col1) && (min(col1)==0 || min(col1)==1) && intlike(col2) && intlike(col3)
            if any(col2 < 0) || any(col3 < 0) || max(abs(col2)) > 5 || max(abs(col3)) > 5
                C.iq = 2; C.jq = 3; C.Re = 6;
                return;
            end
        end
    end
    C.iq = 1; C.jq = 2; C.Re = 5;
end

function print_debug_summary_folderUI_(dbg_out, dbg_cfg, G)
    fprintf("\n================ DEBUG SUMMARY (group14 vs 23) ================\n");
    fprintf("files_total=%d ok=%d path_fail=%d header_fail=%d read_fail=%d cols_fail=%d match_fail=%d | axis=%s\n", ...
        dbg_out.files_total, dbg_out.ok, dbg_out.fail_path, dbg_out.fail_header, dbg_out.fail_read, dbg_out.fail_cols, dbg_out.fail_match, G.u_tag);

    if dbg_cfg.enable
        if ~isempty(dbg_out.bad_path),   fprintf("\n[path_fail examples]\n");   disp(dbg_out.bad_path);   end
        if ~isempty(dbg_out.bad_header), fprintf("\n[header_fail examples]\n"); disp(dbg_out.bad_header); end
        if ~isempty(dbg_out.bad_read),   fprintf("\n[read_fail examples]\n");   disp(dbg_out.bad_read);   end
        if ~isempty(dbg_out.bad_cols),   fprintf("\n[cols_fail examples]\n");   disp(dbg_out.bad_cols);   end
        if ~isempty(dbg_out.bad_match),  fprintf("\n[match_fail examples]\n");  disp(dbg_out.bad_match);  end
    end

    NT = numel(G.T_list); NU = numel(G.U_list);
    if NT==0 || NU==0
        fprintf("[grid] empty.\n");
        fprintf("==================================================================\n\n");
        return;
    end

    miss = ~isfinite(G.chi_map) | ~isfinite(G.doping_map);
    n_miss = nnz(miss);
    fprintf("[grid] missing=%d / %d (%.2f%%)\n", n_miss, NT*NU, 100*n_miss/(NT*NU));

    if n_miss > 0
        [ii,jj] = find(miss);
        fprintf("\n[grid] missing examples as folder pairs:\n");
        for k = 1:min(dbg_cfg.max_missing_print, numel(ii))
            it = ii(k); iu = jj(k);
            ttag = G.Tfolder_map(it,iu);
            utag = G.Ufolder_map(it,iu);
            if ttag==""; ttag = "T"+string(G.T_list(it)); end
            if utag==""; utag = string(G.u_tag)+string(G.U_list(iu)); end
            fprintf("  %s / %s\n", ttag, utag);
        end
    end

    fprintf("==================================================================\n\n");
end
