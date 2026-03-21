function out = phase_diagram_ElectronicPhaseTransition()
% phase_diagram_12vs34_realchi
% Two-color phase diagram:
%   - phase12 (original phase 1 or 2) => red
%   - phase34 (original phase 3 or 4) => blue
%   - invalid/missing => not shown in colorbar
%
% Phase classification:
%   a2crit = (5/6)*b2^2
%   if a2 >= a2crit:
%       a1<0 => phase4 else phase1
%   else:
%       a1<0 => phase3 else phase2
%
% X-axis = folder-controlled variable (mu or doping):
%   - .../Txxx/muXXX/*.txt     => x-axis = mu(folder)
%   - .../Txxx/dopingXXX/*.txt => x-axis = doping(folder)
%
% Internal coef.eval ALWAYS uses doping from HEADER:
%   C = coef.eval(T_header, doping_header, chi_used)

    % -------------------- params / coeff --------------------
    par  = make_realchi_params();
    coef = make_realchi_coeff(par);

    % -------------------- style --------------------
    FS = 18;

    % ---------------- choose folder -----------------
    default_root = "/Users/haoranyan/rg_master/data/";
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

    % ---------------- load grid (folder UI axis + header doping) ----------------
    dbg = struct(); dbg.enable = true; dbg.max_print = 60;
    [G, dbg_out] = load_chi_grid_folderUI_debug_(root, iq_pick, jq_pick, dbg);

    T_list  = G.T_list;      % header T
    U_list  = G.U_list;      % folder mu/doping (display axis)
    u_tag   = G.u_tag;       % 'mu' or 'doping' (axis label)
    chi_map = G.chi_map;     % Re(chi) at (iq,jq)
    dop_map = G.doping_map;  % header doping (for coef.eval)
    Tfolder_map = G.Tfolder_map;
    Ufolder_map = G.Ufolder_map;

    NT = numel(T_list);
    NU = numel(U_list);

    fprintf("Loaded grid: NT=%d, N%s=%d, valid=%d | q=(%d,%d)\n", ...
        NT, u_tag, NU, nnz(isfinite(chi_map) & isfinite(dop_map)), iq_pick, jq_pick);

    if dbg.enable
        print_debug_summary_(dbg_out, G);
    end

    if NT==0 || NU==0
        error("Parsed 0 valid points. Check folder tree / headers / (iq,jq).");
    end

    % ---------------- classify (original 1..4) ----------------
    phase_map = zeros(NT, NU);
    a1_map   = nan(NT, NU);
    a2_map   = nan(NT, NU);
    b2_map   = nan(NT, NU);
    crit_map = nan(NT, NU);

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
            a1_row  = nan(1,NU);
            a2_row  = nan(1,NU);
            b2_row  = nan(1,NU);
            cr_row  = nan(1,NU);

            for iU = 1:NU
                chi_used = chi_map(iT,iU);
                dop = dop_map(iT,iU);

                if ~isfinite(chi_used) || ~isfinite(dop)
                    ph_row(iU) = 0;
                    continue;
                end

                try
                    C  = coef.eval(T, dop, chi_used);
                    a1 = C.a1; a2 = C.a2; b2 = C.b2;
                catch
                    ph_row(iU) = 0;
                    continue;
                end

                a2crit = (5/6) * (b2^2);

                a1_row(iU) = a1;
                a2_row(iU) = a2;
                b2_row(iU) = b2;
                cr_row(iU) = a2crit;

                if ~isfinite(a1) || ~isfinite(a2) || ~isfinite(b2) || ~isfinite(a2crit)
                    ph_row(iU) = 0;
                else
                    if a2 >= a2crit
                        ph_row(iU) = 1;                % phase1
                        if a1 < 0, ph_row(iU) = 4; end % phase4
                    else
                        ph_row(iU) = 2;                % phase2
                        if a1 < 0, ph_row(iU) = 3; end % phase3
                    end
                end
            end

            phase_map(iT,:) = ph_row;
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
                    continue;
                end

                try
                    C  = coef.eval(T, dop, chi_used);
                    a1 = C.a1; a2 = C.a2; b2 = C.b2;
                catch
                    phase_map(iT,iU) = 0;
                    continue;
                end

                a2crit = (5/6) * (b2^2);

                a1_map(iT,iU)   = a1;
                a2_map(iT,iU)   = a2;
                b2_map(iT,iU)   = b2;
                crit_map(iT,iU) = a2crit;

                if ~isfinite(a1) || ~isfinite(a2) || ~isfinite(b2) || ~isfinite(a2crit)
                    phase_map(iT,iU) = 0;
                else
                    if a2 >= a2crit
                        phase_map(iT,iU) = 1;      % phase1
                        if a1 < 0, phase_map(iT,iU) = 4; end
                    else
                        phase_map(iT,iU) = 2;      % phase2
                        if a1 < 0, phase_map(iT,iU) = 3; end
                    end
                end
            end
        end
    end

    % ---------------- compress for display only: phase12 vs phase34 ----------------
    phase2_map = nan(size(phase_map));
    phase2_map(phase_map==1 | phase_map==2) = 1; % phase12
    phase2_map(phase_map==3 | phase_map==4) = 2; % phase34

    % ---------------- plot ----------------
    fig = figure('Color','w','Name','Phase diagram (phase12 vs phase34)');
    ax = axes(fig);

    % colors: only two phases shown in colorbar
    c12 = [0.85 0.20 0.15];       % red
    c34 = [0.20 0.45 0.95];       % blue
    cmap2 = [c12; c34];

    U_edges = centers_to_edges_(U_list);
    T_edges = centers_to_edges_(T_list);

    Cgrid = nan(numel(T_edges), numel(U_edges));
    Cgrid(1:end-1, 1:end-1) = phase2_map;

    [UUe, TTe] = meshgrid(U_edges, T_edges);
    surface(ax, UUe, TTe, zeros(size(UUe)), Cgrid, ...
        'EdgeColor','none', 'FaceColor','flat');
    view(ax, 2);
    axis(ax, 'tight');
    set(ax,'YDir','normal');

    colormap(ax, cmap2);
    caxis(ax, [1 2]);

    ax.XTickMode = 'auto';
    ax.YTickMode = 'auto';

    if strcmpi(u_tag,'doping')
        xlabel(ax, '$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$', ...
            'Interpreter','latex','FontSize',FS);
        xtickformat(ax,'%.4f');
    else
        xlabel(ax, '$\mu\ (\mathrm{eV})$', ...
            'Interpreter','latex','FontSize',FS);
        xtickformat(ax,'%.6g');
    end
    ylabel(ax, '$T$', 'Interpreter','latex','FontSize',FS);

    title(ax, 'Electronic phase transition', ...
        'Interpreter','latex','FontWeight','normal','FontSize',FS);

    set(ax, 'FontSize',FS, ...
        'TickLabelInterpreter','latex', ...
        'LineWidth',1, ...
        'TickDir','out', ...
        'Box','on');

    % colorbar: only red/blue two phases
    cb = colorbar(ax, 'eastoutside');
    cb.Ticks = [1 2];
    cb.TickLabels = {'phase12', 'phase34'};
    cb.FontSize = FS;
    cb.TickLabelInterpreter = 'none';
    cb.Label.String = 'Phase';
    cb.Label.FontSize = FS;
    cb.Label.Interpreter = 'none';

    % datatip
    dcm = datacursormode(fig);
    set(dcm, 'Enable','on', 'SnapToDataVertex','off', 'DisplayStyle','datatip');
    set(dcm, 'UpdateFcn', @tip_cb_);

    % output struct
    out = struct();
    out.phase_map   = phase_map;    % original 0..4
    out.phase2_map  = phase2_map;   % display 1..2 / NaN
    out.T_list      = T_list;
    out.U_list      = U_list;
    out.u_tag       = u_tag;
    out.chi_map     = chi_map;
    out.doping_map  = dop_map;
    out.a1_map      = a1_map;
    out.a2_map      = a2_map;
    out.b2_map      = b2_map;
    out.crit_map    = crit_map;
    out.Tfolder_map = Tfolder_map;
    out.Ufolder_map = Ufolder_map;
    out.iq_pick     = iq_pick;
    out.jq_pick     = jq_pick;

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

        ph0  = phase_map(it, iu);
        grp  = phase2_map(it, iu);
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

        gname = "invalid";
        if isequal(grp,1), gname="phase12"; end
        if isequal(grp,2), gname="phase34"; end

        txt = {
            sprintf('T = %.6g K', T0)
            Ustr
            sprintf('group = %s', gname)
            sprintf('phase(original) = %d (%s)', ph0, phase_name1234_(ph0))
            sprintf('doping(header) = %.6g', dop0)
            sprintf('a1 = %.6g', a10)
            sprintf('a2 = %.6g', a20)
            sprintf('b2 = %.6g', b20)
            sprintf('a2crit = %.6g', cr0)
            sprintf('folder = %s / %s', string(Tfolder_map(it,iu)), string(Ufolder_map(it,iu)))
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

% =====================================================================
% Loader: x-axis from folder (mu/doping), T from header, doping from header, chi from table
% =====================================================================
function [G, dbg] = load_chi_grid_folderUI_debug_(root_dir, iq_pick, jq_pick, dbg_cfg)
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

    Tfolder = strings(n,1);
    Ufolder = strings(n,1);

    u_tag = "";

    n_ok=0; n_fail_path=0; n_fail_header=0; n_fail_read=0; n_fail_cols=0; n_fail_match=0;
    bad_path=strings(0); bad_header=strings(0); bad_read=strings(0); bad_cols=strings(0); bad_match=strings(0);

    for k = 1:n
        fpath = string(fullfile(L(k).folder, L(k).name));

        [Tfd, Ufd, ukind, uval] = parse_T_and_U_from_path_(fpath);
        Tfolder(k) = Tfd;
        Ufolder(k) = Ufd;

        if u_tag == "" && ukind ~= ""
            u_tag = ukind;
        end

        if ukind == "" || ~isfinite(uval)
            n_fail_path = n_fail_path + 1;
            if dbg_cfg.enable && numel(bad_path) < dbg_cfg.max_print
                bad_path(end+1) = sprintf("%s | cannot parse folder U | %s/%s", fpath, Tfd, Ufd);
            end
            continue;
        end

        H = parse_header_T_and_doping_robust_(fpath);
        if ~H.ok
            n_fail_header = n_fail_header + 1;
            if dbg_cfg.enable && numel(bad_header) < dbg_cfg.max_print
                bad_header(end+1) = sprintf("%s | header missing T/doping (T=%g doping=%g) | %s/%s", ...
                    fpath, H.T_K, H.doping, Tfd, Ufd);
            end
            continue;
        end

        try
            M = read_chi_numeric_skiphash_(fpath);
        catch ME
            n_fail_read = n_fail_read + 1;
            if dbg_cfg.enable && numel(bad_read) < dbg_cfg.max_print
                bad_read(end+1) = sprintf("%s | read failed: %s | %s/%s", fpath, string(ME.message), Tfd, Ufd);
            end
            continue;
        end

        if isempty(M) || size(M,2) < 6
            n_fail_cols = n_fail_cols + 1;
            if dbg_cfg.enable && numel(bad_cols) < dbg_cfg.max_print
                bad_cols(end+1) = sprintf("%s | numeric table too small: %dx%d | %s/%s", ...
                    fpath, size(M,1), size(M,2), Tfd, Ufd);
            end
            continue;
        end

        Cc = detect_cols_(M);
        id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
        if isempty(id)
            n_fail_match = n_fail_match + 1;
            if dbg_cfg.enable && numel(bad_match) < dbg_cfg.max_print
                bad_match(end+1) = sprintf("%s | no (iq,jq)=(%d,%d) | %s/%s", fpath, iq_pick, jq_pick, Tfd, Ufd);
            end
            continue;
        end

        ReChi = M(id, Cc.Re);

        Tvals(k) = round(H.T_K, 12);
        Uvals(k) = quantize_U_(ukind, uval);
        chiV(k)  = ReChi;
        doping_header(k) = round(H.doping, 12);

        n_ok = n_ok + 1;
    end

    if u_tag == ""
        u_tag = "doping";
    end

    dbg = struct();
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
    Tfolder = Tfolder(mask);
    Ufolder = Ufolder(mask);

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
            Tfolder_map(iT,iU) = Tfolder(i);
            Ufolder_map(iT,iU) = Ufolder(i);
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

    extract_num = @(s) str2double(regexp(char(s), '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)', 'match', 'once'));

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

function print_debug_summary_(dbg, G)
    fprintf("\n================ DEBUG SUMMARY (load folderUI) ================\n");
    fprintf("files_total=%d ok=%d path_fail=%d header_fail=%d read_fail=%d cols_fail=%d match_fail=%d | axis=%s\n", ...
        dbg.files_total, dbg.ok, dbg.fail_path, dbg.fail_header, dbg.fail_read, dbg.fail_cols, dbg.fail_match, G.u_tag);

    NT = numel(G.T_list); NU = numel(G.U_list);
    if NT==0 || NU==0
        fprintf("[grid] empty.\n");
        fprintf("================================================================\n\n");
        return;
    end
    miss = ~isfinite(G.chi_map) | ~isfinite(G.doping_map);
    n_miss = nnz(miss);
    fprintf("[grid] missing=%d / %d (%.2f%%)\n", n_miss, NT*NU, 100*n_miss/(NT*NU));

    fprintf("================================================================\n\n");
end