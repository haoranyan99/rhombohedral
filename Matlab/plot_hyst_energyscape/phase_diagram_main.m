function out = phase_diagram_main()
% phase_diagram_realchi4_folderUI_aligned_hover (SIMPLIFIED FIX)
% Fixes:
%   1) remove dependency on par.fontSize (use FS)
%   2) title uses Interpreter='none' to avoid LaTeX syntax errors
%   3) delete helper phase_name_ (datatip shows phase number only)

    % -------------------- params / coeff --------------------
    par  = make_realchi_params();
    coef = make_realchi_coeff(par);

    % ---------------- choose folder -----------------
    default_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_200\D0.084";
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

    % ---------------- classification settings (YOUR v3 rules) ----------------
    cls = struct();
    cls.h = 2e-3;                 % FD step for Hessian
    cls.eig_eps = 1e-10;          % PD threshold
    cls.X0_tol = 1;               % |X|>X0_tol treated as X != 0
    cls.cluster_tol_1d = 3e-2;    % merge nearby 1D minima locations

    cls.origin_box = 1;
    cls.origin_grid = 11;
    cls.origin_refine_topK = 6;
    cls.origin_refine_jitter = 3e-5;

    cls.Npsi_small = 7;
    cls.Nx_small   = 401;

    cls.Npsi_wide = 41;
    cls.Nx_wide   = 601;
    cls.psi_line_max = 0.9 * max(abs(par.psi1_lim));

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

    % parallel (optional)
    try
        p = gcp('nocreate');
        if isempty(p), parpool; end
    catch
        fprintf("[warn] Parallel pool not available. Running serial.\n");
    end

    % ---------------- classify ----------------
    phase_map = zeros(NT, NU);

    parfor iT = 1:NT
        T = T_list(iT);
        ph_row = zeros(1,NU);

        for iU = 1:NU
            chi_used = chi_map(iT,iU);
            dop = dop_map(iT,iU);  % header doping used internally

            if ~isfinite(chi_used) || ~isfinite(dop)
                ph_row(iU) = 0;
                continue;
            end

            C = coef.eval(T, dop, chi_used);

            try
                F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
            catch
                ph_row(iU) = 0;
                continue;
            end

            ph_row(iU) = classify_point_(F, par.psi1_lim, par.psi2_lim, cls);
        end

        phase_map(iT,:) = ph_row;
    end

    % ---------------- plot (aligned blocks + datatip) ----------------
    FS = 14; % <--- simplified: fixed font size, no par.fontSize
    fig = figure('Color','w','Name','Phase diagram (folder-UI axis, aligned)');
    ax = axes(fig); %#ok<LAXES>

    % colors: [0..4]
    c_invalid = [0.25 0.08 0.35];
    c1 = [0.20 0.55 0.95];
    c2 = [0.25 0.85 0.35];
    c3 = [1.00 0.60 0.10];
    c4 = [0.85 0.20 0.15];
    cmap = [c_invalid; c1; c2; c3; c4];

    % compute edges so each block aligns with true coordinate bins
    U_edges = centers_to_edges_(U_list);
    T_edges = centers_to_edges_(T_list);

    % C for surface: needs (NT+1) x (NU+1)
    Cgrid = nan(numel(T_edges), numel(U_edges));
    Cgrid(1:end-1, 1:end-1) = phase_map;

    [UU, TT] = meshgrid(U_edges, T_edges);

    surface(ax, UU, TT, zeros(size(UU)), Cgrid, ...
        'EdgeColor','none', 'FaceColor','flat');
    view(ax, 2);
    axis(ax, 'tight');
    set(ax,'YDir','normal');

    colormap(ax, cmap);
    caxis(ax, [0 4]);

    % ticks at centers
    xticks(ax, U_list);
    yticks(ax, T_list);

    if strcmpi(u_tag,'doping')
        xlabel(ax, '$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$','Interpreter','latex');
    else
        xlabel(ax, '$\mu\ (\mathrm{eV})$','Interpreter','latex');
    end
    ylabel(ax, '$T$','Interpreter','latex');

    % <--- simplified: title uses Interpreter='none' (no LaTeX parse errors)
    title(ax, sprintf('Phase map (real chi), q=(%d,%d), lambda=%.3g | x-axis=%s(folder)', ...
        iq_pick, jq_pick, par.lambda, u_tag), ...
        'Interpreter','none','FontWeight','normal');

    set(ax,'FontSize',FS,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out','Box','on');

    % legend (no colorbar)
    hold(ax,'on');
    h0 = plot(ax, nan,nan,'s','MarkerFaceColor',c_invalid,'MarkerEdgeColor','none','MarkerSize',10);
    h1 = plot(ax, nan,nan,'s','MarkerFaceColor',c1,'MarkerEdgeColor','none','MarkerSize',10);
    h2 = plot(ax, nan,nan,'s','MarkerFaceColor',c2,'MarkerEdgeColor','none','MarkerSize',10);
    h3 = plot(ax, nan,nan,'s','MarkerFaceColor',c3,'MarkerEdgeColor','none','MarkerSize',10);
    h4 = plot(ax, nan,nan,'s','MarkerFaceColor',c4,'MarkerEdgeColor','none','MarkerSize',10);
    hold(ax,'off');

    legend([h1 h2 h3 h4 h0], ...
        {'phase1','phase2','phase3','phase4','invalid / missing'}, ...
         'Interpreter','tex','Location','eastoutside');

    % datatip / hover
    dcm = datacursormode(fig);
    set(dcm, 'Enable','on', 'SnapToDataVertex','off', 'DisplayStyle','datatip');
    set(dcm, 'UpdateFcn', @tip_cb_);

    % ---------------- save ----------------
    out_dir = fullfile(root, "..", "plot");
    if ~exist(out_dir,"dir"), mkdir(out_dir); end
    out_png = fullfile(out_dir, sprintf("phase_realchi_folderUI_v4_aligned_iq%d_jq%d_lambda%.3g_x%s.png", ...
        iq_pick, jq_pick, par.lambda, u_tag));
    exportgraphics(fig, out_png, 'Resolution', 300);
    fprintf("Saved phase: %s\n", out_png);

    % output struct
    out = struct();
    out.phase_map = phase_map;
    out.T_list = T_list;
    out.U_list = U_list;
    out.u_tag = u_tag;
    out.chi_map = chi_map;
    out.doping_map = dop_map;
    out.iq_pick = iq_pick;
    out.jq_pick = jq_pick;
    out.lambda = par.lambda;
    out.cls = cls;
    out.png = out_png;

    % ===================== nested: datatip callback =====================
    function txt = tip_cb_(~, evt)
        pos = evt.Position;   % [x y z]
        x = pos(1);
        y = pos(2);

        iu = find(x >= U_edges(1:end-1) & x < U_edges(2:end), 1, 'first');
        it = find(y >= T_edges(1:end-1) & y < T_edges(2:end), 1, 'first');

        if isempty(iu) || isempty(it)
            txt = {'(out of range)'};
            return;
        end

        ph = phase_map(it, iu);

        T0 = T_list(it);
        U0 = U_list(iu);
        dop0 = dop_map(it, iu);

        if strcmpi(u_tag,'doping')
            Ustr = sprintf('doping(folder) = %.4f', U0);
        else
            Ustr = sprintf('mu(folder) = %.6g', U0);
        end

        txt = {
            sprintf('T = %.6g K', T0)
            Ustr
            sprintf('phase = %d', ph)          % <--- simplified: no phase_name_ helper
            sprintf('doping = %.6g', dop0)
        };
    end
end

% =====================================================================
% One-point classifier (YOUR final rules, unchanged)
% =====================================================================
function ph = classify_point_(F, psi_lim, X_lim, cls)
    ph = 0;

    try
        v0 = F(0,0);
        if ~isfinite(v0), return; end
    catch
        return;
    end

    ok0 = exists_local_min_in_box_(F, cls.origin_box, cls.origin_grid, cls.origin_refine_topK, ...
                                  cls.origin_refine_jitter, cls.h, cls.eig_eps, psi_lim, X_lim);

    if ok0
        psi_list = linspace(-cls.origin_box, cls.origin_box, cls.Npsi_small);
        hasX = false;
        for k = 1:numel(psi_list)
            psi0 = psi_list(k);
            if psi0 < psi_lim(1) || psi0 > psi_lim(2), continue; end
            if exists_nonzero_min_on_line_( @(X) F(psi0,X), X_lim, cls.Nx_small, cls.cluster_tol_1d, cls.X0_tol)
                hasX = true;
                break;
            end
        end
        ph = 1 + hasX;
        return;
    end

    psi_list = linspace(-cls.psi_line_max, cls.psi_line_max, cls.Npsi_wide);
    for k = 1:numel(psi_list)
        psi0 = psi_list(k);
        if psi0 < psi_lim(1) || psi0 > psi_lim(2), continue; end
        if exists_nonzero_min_on_line_( @(X) F(psi0,X), X_lim, cls.Nx_wide, cls.cluster_tol_1d, cls.X0_tol)
            ph = 4;
            return;
        end
    end

    ph = 3;
end

function ok = exists_local_min_in_box_(F, box, Ng, topK, jitter, h, eig_eps, psi_lim, X_lim)
    ok = false;

    box_psi = [max(-box, psi_lim(1)), min(+box, psi_lim(2))];
    box_X   = [max(-box, X_lim(1)),   min(+box, X_lim(2))];

    psi = linspace(box_psi(1), box_psi(2), Ng);
    X   = linspace(box_X(1),   box_X(2),   Ng);
    [P, Q] = meshgrid(psi, X);

    V = arrayfun(@(u,v) safe_eval_(F,u,v), P, Q);
    V(~isfinite(V)) = 1e30;

    mask = local_min_mask_2d_(V);
    [rr,cc] = find(mask);
    if isempty(rr), return; end

    cand = [P(sub2ind(size(P), rr, cc)), Q(sub2ind(size(Q), rr, cc))];
    Ec   = V(sub2ind(size(V), rr, cc));
    [~,ix] = sort(Ec,'ascend');
    ix = ix(1:min(topK, numel(ix)));
    cand = cand(ix,:);

    fopt = optimset('Display','off','MaxIter',200,'TolX',1e-10,'TolFun',1e-12);
    for i = 1:size(cand,1)
        x0 = cand(i,:) + jitter*randn(1,2);
        x0(1) = min(max(x0(1), box_psi(1)), box_psi(2));
        x0(2) = min(max(x0(2), box_X(1)),   box_X(2));

        try
            xhat = fminsearch(@(x) safe_eval_(F,x(1),x(2)), x0, fopt);
        catch
            continue;
        end

        if ~isfinite(xhat(1)) || ~isfinite(xhat(2)), continue; end
        if xhat(1) < box_psi(1) || xhat(1) > box_psi(2), continue; end
        if xhat(2) < box_X(1)   || xhat(2) > box_X(2),   continue; end

        H = hessian_fd_(@(u,v) F(u,v), xhat(1), xhat(2), h);
        ev = eig((H+H')/2);
        if all(real(ev) > eig_eps)
            ok = true;
            return;
        end
    end
end

function hasNonzero = exists_nonzero_min_on_line_(f1d, xlim, N, merge_tol, X0_tol)
    hasNonzero = false;

    x = linspace(xlim(1), xlim(2), N);
    y = arrayfun(@(t) safe_eval_1d_(f1d,t), x);
    y(~isfinite(y)) = 1e30;

    if N < 3, return; end
    ym1 = y(1:end-2);
    y0  = y(2:end-1);
    yp1 = y(3:end);
    mask = (y0 <= ym1) & (y0 <= yp1);
    idx = find(mask) + 1;
    if isempty(idx), return; end

    xc = sort(x(idx));
    xuniq = xc(1);
    for i = 2:numel(xc)
        if abs(xc(i) - xuniq(end)) > merge_tol
            xuniq(end+1) = xc(i); %#ok<AGROW>
        end
    end

    hasNonzero = any(abs(xuniq) > X0_tol);
end

function v = safe_eval_(F, psi, X)
    try
        v = F(psi,X);
        if ~isfinite(v), v = 1e30; end
    catch
        v = 1e30;
    end
end

function v = safe_eval_1d_(f, x)
    try
        v = f(x);
        if ~isfinite(v), v = 1e30; end
    catch
        v = 1e30;
    end
end

function mask = local_min_mask_2d_(A)
    [R,C] = size(A);
    mask = false(R,C);
    if R < 3 || C < 3, return; end

    center = A(2:R-1, 2:C-1);
    n1 = A(1:R-2, 2:C-1);
    n2 = A(3:R  , 2:C-1);
    n3 = A(2:R-1, 1:C-2);
    n4 = A(2:R-1, 3:C  );
    n5 = A(1:R-2, 1:C-2);
    n6 = A(1:R-2, 3:C  );
    n7 = A(3:R  , 1:C-2);
    n8 = A(3:R  , 3:C  );

    m = center <= n1 & center <= n2 & center <= n3 & center <= n4 & ...
        center <= n5 & center <= n6 & center <= n7 & center <= n8;

    mask(2:R-1, 2:C-1) = m;
end

function H = hessian_fd_(f, x, y, h)
    f00 = f(x,y);
    fxx = (f(x+h,y) - 2*f00 + f(x-h,y))/h^2;
    fyy = (f(x,y+h) - 2*f00 + f(x,y-h))/h^2;
    fxy = (f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)) / (4*h^2);
    H = [fxx, fxy; fxy, fyy];
end

% =====================================================================
% Folder-UI loader (mu/doping from path), doping from header, chi from table
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
    fprintf("\n================ DEBUG SUMMARY (phase v4 folderUI) ================\n");
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

% =====================================================================
% helpers: edges
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