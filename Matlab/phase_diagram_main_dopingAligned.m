function out = phase_diagram_main_dopingAligned()
% Phase diagram with X-axis = doping (from HEADER), aligned across T.
% Folder mu/doping tag is optional only for datatip display.

    par  = make_realchi_params();
    coef = make_realchi_coeff(par);

    default_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_200\D0.067";
    if ~isfolder(default_root), default_root = string(pwd); end
    root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
    if isequal(root,0), error("User cancelled."); end
    root = string(root);

    iq_pick = round(par.iq_pick);
    jq_pick = round(par.jq_pick);

    % ---------- load points ----------
    P = load_chi_points_headerDoping_(root, iq_pick, jq_pick);

    if isempty(P.T)
        error("No valid chi points loaded. Check headers / (iq,jq) / files.");
    end

    % ---------- classify each point ----------
    cls = make_cls_defaults_(par); % 你也可以把这块替换成你现在脚本里的 cls
    N = numel(P.T);
    phase_pt = zeros(N,1);

    % optional parallel
    try
        p = gcp('nocreate');
        if isempty(p), parpool; end
    catch
    end

    parfor k = 1:N
        T   = P.T(k);
        dop = P.doping(k);
        chi = P.chi(k);

        if ~isfinite(T) || ~isfinite(dop) || ~isfinite(chi)
            phase_pt(k) = 0; %#ok<PFBNS>
            continue;
        end

        C = coef.eval(T, dop, chi);
        try
            F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
            phase_pt(k) = classify_point_(F, par.psi1_lim, par.psi2_lim, cls);
        catch
            phase_pt(k) = 0;
        end
    end

    % ---------- build aligned grid (T x doping_grid) ----------
    T_list = sort(unique(round(P.T,12)));
    % 全局 doping 网格：用所有点的 unique（更忠实），也可以换成 linspace
    dop_grid = sort(unique(round(P.doping,12)));

    NT = numel(T_list);
    ND = numel(dop_grid);

    phase_map = nan(NT, ND);
    mu_map    = nan(NT, ND);   % 仅用于 datatip（最近邻对应的 mu_folder）
    chi_map   = nan(NT, ND);   % 可选

    for iT = 1:NT
        T0 = T_list(iT);

        idx = find(abs(P.T - T0) < 1e-9 & isfinite(P.doping) & isfinite(phase_pt));
        if isempty(idx), continue; end

        dop_i = P.doping(idx);
        ph_i  = phase_pt(idx);
        mu_i  = P.mu_folder(idx);
        chi_i = P.chi(idx);

        % 对齐到全局 dop_grid：最近邻（相是离散的，推荐 nearest）
        for j = 1:ND
            [dmin, ii] = min(abs(dop_i - dop_grid(j)));
            if ~isfinite(dmin), continue; end
            % 你可以设一个容忍阈值，避免把很远的点硬贴上去
            if dmin > 1e-3    % <- 这里你按自己 doping 分辨率改
                continue;
            end
            phase_map(iT,j) = ph_i(ii);
            mu_map(iT,j)    = mu_i(ii);
            chi_map(iT,j)   = chi_i(ii);
        end
    end

    % ---------- plot aligned blocks ----------
    fig = figure('Color','w','Name','Phase diagram (aligned doping)');
    ax = axes(fig);

    dop_edges = centers_to_edges_(dop_grid);
    T_edges   = centers_to_edges_(T_list);

    Cgrid = nan(numel(T_edges), numel(dop_edges));
    Cgrid(1:end-1,1:end-1) = phase_map;

    [XX, YY] = meshgrid(dop_edges, T_edges);

    surface(ax, XX, YY, zeros(size(XX)), Cgrid, 'EdgeColor','none','FaceColor','flat');
    view(ax,2); axis(ax,'tight'); set(ax,'YDir','normal');

    cmap = [0.25 0.08 0.35; 0.20 0.55 0.95; 0.25 0.85 0.35; 1.00 0.60 0.10; 0.85 0.20 0.15];
    colormap(ax, cmap); caxis(ax,[0 4]);

    xlabel(ax,'$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$','Interpreter','latex');
    ylabel(ax,'$T$','Interpreter','latex');
    title(ax, sprintf('Phase map (real \\chi), q=(%d,%d), \\lambda=%.3g | x=doping(header)', ...
        iq_pick, jq_pick, par.lambda), 'Interpreter','latex','FontWeight','normal');

    set(ax,'FontSize',par.fontSize,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out','Box','on');

    % datatip: 显示 doping + 对应的 mu(folder)
    dcm = datacursormode(fig);
    set(dcm,'Enable','on','SnapToDataVertex','off','DisplayStyle','datatip');
    set(dcm,'UpdateFcn', @tip_cb_);

    out_dir = fullfile(root, "..", "plot");
    if ~exist(out_dir,"dir"), mkdir(out_dir); end
    out_png = fullfile(out_dir, sprintf("phase_realchi_dopingAligned_iq%d_jq%d_lambda%.3g.png", ...
        iq_pick, jq_pick, par.lambda));
    exportgraphics(fig, out_png, 'Resolution', 300);

    out = struct();
    out.T_list = T_list;
    out.dop_grid = dop_grid;
    out.phase_map = phase_map;
    out.mu_map = mu_map;
    out.chi_map = chi_map;
    out.png = out_png;

    % ---------- nested datatip ----------
    function txt = tip_cb_(~, evt)
        pos = evt.Position;
        x = pos(1); y = pos(2);
        j = find(x >= dop_edges(1:end-1) & x < dop_edges(2:end), 1, 'first');
        i = find(y >= T_edges(1:end-1) & y < T_edges(2:end), 1, 'first');
        if isempty(i) || isempty(j)
            txt = {'(out of range)'}; return;
        end
        ph = phase_map(i,j);
        txt = {
            sprintf('T = %.6g K', T_list(i))
            sprintf('doping = %.6g', dop_grid(j))
            sprintf('phase = %g', ph)
            sprintf('mu(folder, nearest) = %.6g', mu_map(i,j))
            sprintf('Re chi = %.6g', chi_map(i,j))
        };
    end
end

% ===== load points: (T, doping(header), chi(iq,jq), mu(folder optional)) =====
function P = load_chi_points_headerDoping_(root_dir, iq_pick, jq_pick)
    L = dir(fullfile(root_dir, "**", "chi*.txt"));
    if isempty(L)
        P = struct('T',[],'doping',[],'chi',[],'mu_folder',[]);
        return;
    end

    T = []; doping = []; chi = []; mu_folder = [];

    for k = 1:numel(L)
        fpath = fullfile(L(k).folder, L(k).name);

        H = parse_header_T_and_doping_robust_(fpath);
        if ~H.ok, continue; end

        try
            M = read_chi_numeric_skiphash_(fpath);
        catch
            continue;
        end
        if isempty(M) || size(M,2) < 6, continue; end

        Cc = detect_cols_(M);
        id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
        if isempty(id), continue; end

        ReChi = M(id, Cc.Re);
        T(end+1,1) = H.T_K; %#ok<AGROW>
        doping(end+1,1) = H.doping; %#ok<AGROW>
        chi(end+1,1) = ReChi; %#ok<AGROW>

        % optional: parse mu folder for info only
        mu_folder(end+1,1) = parse_mu_from_path_(string(fpath)); %#ok<AGROW>
    end

    P = struct('T',T,'doping',doping,'chi',chi,'mu_folder',mu_folder);
end

function mu = parse_mu_from_path_(fpath)
    mu = NaN;
    p = replace(string(fpath),"\","/");
    tok = regexp(p,'/(mu[^/]+)/','tokens','once');
    if isempty(tok), return; end
    s = string(tok{1});
    num = regexp(s,'([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)','match','once');
    if ~isempty(num), mu = str2double(num); end
end

function cls = make_cls_defaults_(par)
    cls = struct();
    cls.h = 2e-3;
    cls.eig_eps = 1e-10;
    cls.X0_tol = 1;
    cls.cluster_tol_1d = 3e-2;
    cls.origin_box = 1;
    cls.origin_grid = 11;
    cls.origin_refine_topK = 6;
    cls.origin_refine_jitter = 3e-5;
    cls.Npsi_small = 7;
    cls.Nx_small = 401;
    cls.Npsi_wide = 41;
    cls.Nx_wide = 601;
    cls.psi_line_max = 0.9 * max(abs(par.psi1_lim));
end

% ===== you already have these in your file; keep or reuse =====
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
    fprintf("\n================ DEBUG SUMMARY (phase a2crit + a1 split) ================\n");
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