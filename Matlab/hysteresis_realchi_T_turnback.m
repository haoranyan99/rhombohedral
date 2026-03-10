function H = hysteresis_realchi_T_turnback()
% hysteresis_realchi_T_turnback
% Hysteresis at ONE temperature: sweep doping forward to a TURN point,
% then immediately sweep backward (return scan).
%
% Unified convention:
%   - absolute (iq_pick, jq_pick)
%   - chi_used = Re(chi) (no sign flip)
%   - C = coef.eval(T, dop, chi_used) returns a1,b1,a2,b2,c2
%
% Return-scan logic:
%   Stage-1: dop: dmin -> d_turn (N1 points)
%   Stage-2: dop: d_turn -> dmin (N2 points)
%
% Continuation:
%   pick minimum closest to previous solution.

    % ---------------- params / coeff ----------------
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

    % ---------------- q point (ABS) -----------------
    iq_pick = round(par.iq_pick);
    jq_pick = round(par.jq_pick);
    fprintf("Pick absolute q-point: (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);

    % ---------------- load chi grid -----------------
    D = load_chi_grid_newpatch_absq_recursive_(root, iq_pick, jq_pick);
    T_list   = D.T_list;
    dop_list = D.dop_list;
    chi_map  = D.chi_used_map;

    NT = numel(T_list);
    Nd = numel(dop_list);
    fprintf("Loaded grid: NT=%d, Ndop=%d, valid=%d\n", NT, Nd, nnz(isfinite(chi_map)));

    if NT==0 || Nd==0
        error("Parsed 0 valid points.");
    end

    % ---------------- choose T ----------------------
    if isfield(par,'T_target')
        T_target = par.T_target;
    else
        T_target = T_list(1);
    end
    [~, iT0] = min(abs(T_list - T_target));
    T0 = T_list(iT0);
    fprintf("Return-scan hysteresis at T_target=%.12g -> using T=%.12g (index %d/%d)\n", ...
        T_target, T0, iT0, NT);

    % valid dop points at this T
    chi_row = chi_map(iT0,:).';
    mask = isfinite(chi_row);
    dop0 = dop_list(:);
    dop0 = dop0(mask);
    chi0 = chi_row(mask);

    if numel(dop0) < 3
        error("Too few valid dop points at this T to do hysteresis.");
    end

    % sort by doping
    [dop0, ord] = sort(dop0);
    chi0 = chi0(ord);

    % chi interpolation function
    chi_of_dop = @(d) interp1(dop0, chi0, d, 'linear', 'extrap');

    % ---------------- return-scan grid --------------
    dmin = min(dop0);
    dmax = max(dop0);

    % user-specified turn point
    if isfield(par,'hyst') && isfield(par.hyst,'turn_dop') && isfinite(par.hyst.turn_dop)
        d_turn = par.hyst.turn_dop;
    else
        d_turn = 0.5*(dmin + dmax); % default: mid point
    end
    d_turn = min(max(d_turn, dmin), dmax);

    % total points (reuse par.hyst.N if exists)
    if isfield(par,'hyst') && isfield(par.hyst,'N') && isfinite(par.hyst.N) && par.hyst.N>=10
        Ntot = round(par.hyst.N);
    else
        Ntot = 200;
    end

    % split points between forward/backward (include turn point once)
    N1 = max(5, round(0.5*Ntot));
    N2 = Ntot - N1 + 1;

    dop_fwd = linspace(dmin, d_turn, N1).';
    dop_bwd = linspace(d_turn, dmin, N2).';

    fprintf("Return scan: dmin=%.6g -> d_turn=%.6g -> dmin, N1=%d N2=%d\n", dmin, d_turn, N1, N2);

    % ---------------- minima-search settings ---------
    scan = struct();
    scan.Npsi1 = 81;
    scan.Npsi2 = 81;
    scan.use_gaussian_smooth = true;
    scan.smooth_sigma = 0.7;
    scan.keep_topK = 40;

    opt = struct();
    opt.seed_jitter  = 0.12;
    opt.seed_repeats = 2;
    opt.max_iter     = 450;
    opt.tol_fun      = 1e-10;
    opt.tol_x        = 1e-8;

    opt.cluster_tol  = 3e-2;
    opt.fd_h         = 2e-3;
    opt.min_eig_eps  = 1e-6;

    opt.force_origin_and_nearby = true;
    opt.nearby_delta = 0.2;

    % ---------------- stage-1 forward ----------------
    psi_f = nan(N1,2);
    a1_f  = nan(N1,1);
    chi_f = nan(N1,1);

    if isfield(par,'hyst') && isfield(par.hyst,'psi1_0') && isfield(par.hyst,'psi2_0')
        prev = [par.hyst.psi1_0, par.hyst.psi2_0];
    else
        prev = [0,0];
    end

    for i = 1:N1
        dop = dop_fwd(i);
        chi_used = chi_of_dop(dop);
        C = coef.eval(T0, dop, chi_used);

        F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
        mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);

        if isempty(mins)
            [p1,p2,~] = minimize_free_energy(C.a1,C.b1,C.a2,C.b2,C.c2,par.lambda, prev(1),prev(2), false);
            pick = [p1,p2];
        else
            pick = pick_closest_min_(mins, prev);
        end

        psi_f(i,:) = pick;
        prev = pick;

        a1_f(i)  = C.a1;
        chi_f(i) = chi_used;
    end

    % ---------------- stage-2 backward (return) ------
    psi_b = nan(N2,2);
    a1_b  = nan(N2,1);
    chi_b = nan(N2,1);

    prev = psi_f(end,:); % start from turn point state

    for i = 1:N2
        dop = dop_bwd(i);
        chi_used = chi_of_dop(dop);
        C = coef.eval(T0, dop, chi_used);

        F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
        mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);

        if isempty(mins)
            [p1,p2,~] = minimize_free_energy(C.a1,C.b1,C.a2,C.b2,C.c2,par.lambda, prev(1),prev(2), false);
            pick = [p1,p2];
        else
            pick = pick_closest_min_(mins, prev);
        end

        psi_b(i,:) = pick;
        prev = pick;

        a1_b(i)  = C.a1;
        chi_b(i) = chi_used;
    end

    % ---------------- plot --------------------------
    fig = figure('Color','w','Units','pixels','Position',[80 80 1200 720], ...
        'Name','Hysteresis (return scan, real chi, dense minima)');

    tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    fs = par.fontSize;

    % chi(dop)
    ax1 = nexttile(tl,1);
    plot(ax1, dop0, chi0, 'o-','LineWidth',2,'MarkerSize',6);
    grid(ax1,'on'); box(ax1,'on');
    xlabel(ax1,'$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$','Interpreter','latex');
    ylabel(ax1,'$\chi_{\mathrm{used}}=\mathrm{Re}\,\chi$','Interpreter','latex');
    title(ax1, sprintf('$T=%.6g$, q=(%d,%d), turn=%.6g', T0, iq_pick, jq_pick, d_turn), ...
        'Interpreter','latex','FontWeight','normal','FontSize',fs);
    set(ax1,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out');
    xline(ax1, d_turn, '--', 'LineWidth', 1.3);

    % a1(dop) forward/back
    ax2 = nexttile(tl,3);
    plot(ax2, dop_fwd, a1_f, 'k-','LineWidth',2); hold(ax2,'on');
    plot(ax2, dop_bwd, a1_b, 'k--','LineWidth',2); hold(ax2,'off');
    grid(ax2,'on'); box(ax2,'on');
    yline(ax2,0,'--','LineWidth',1.2);
    xlabel(ax2,'$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$','Interpreter','latex');
    ylabel(ax2,'$a_1$','Interpreter','latex');
    title(ax2, sprintf('$a_1=s_\\chi(\\mathrm{invV}-\\chi_{used})$, $s_\\chi=%.3g$, invV=%.6g', ...
        par.chi_scaler, par.invV), ...
        'Interpreter','latex','FontWeight','normal','FontSize',fs);
    legend(ax2,{'forward','return'},'Interpreter','latex','Location','best');
    set(ax2,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out');

    % |psi|
    ax3 = nexttile(tl,2);
    plot(ax3, dop_fwd, abs(psi_f(:,1)), '-','LineWidth',2); hold(ax3,'on');
    plot(ax3, dop_bwd, abs(psi_b(:,1)), '--','LineWidth',2); hold(ax3,'off');
    grid(ax3,'on'); box(ax3,'on');
    xlabel(ax3,'$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$','Interpreter','latex');
    ylabel(ax3,'$|\psi|$','Interpreter','latex');
    title(ax3,'$|\psi|$ (forward then return)','Interpreter','latex','FontWeight','normal','FontSize',fs);
    legend(ax3,{'forward','return'},'Interpreter','latex','Location','best');
    set(ax3,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out');

    % |X|
    ax4 = nexttile(tl,4);
    plot(ax4, dop_fwd, abs(psi_f(:,2)), '-','LineWidth',2); hold(ax4,'on');
    plot(ax4, dop_bwd, abs(psi_b(:,2)), '--','LineWidth',2); hold(ax4,'off');
    grid(ax4,'on'); box(ax4,'on');
    xlabel(ax4,'$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$','Interpreter','latex');
    ylabel(ax4,'$|X|$','Interpreter','latex');
    title(ax4,'$|X|$ (forward then return)','Interpreter','latex','FontWeight','normal','FontSize',fs);
    legend(ax4,{'forward','return'},'Interpreter','latex','Location','best');
    set(ax4,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out');

    linkaxes([ax1 ax2 ax3 ax4],'x');

    out_dir = fullfile(root, "..", "plot");
    if ~exist(out_dir,"dir"), mkdir(out_dir); end
    out_png = fullfile(out_dir, sprintf("hyst_return_realchi_T%.6g_turn%.6g_iq%d_jq%d_lambda%.3g.png", ...
        T0, d_turn, iq_pick, jq_pick, par.lambda));
    exportgraphics(fig, out_png, 'Resolution', 300);
    fprintf("Saved return-scan hyst: %s\n", out_png);

    % ---------------- output struct ----------------
    H = struct();
    H.root = root;
    H.iq_pick = iq_pick;
    H.jq_pick = jq_pick;
    H.T = T0;
    H.lambda = par.lambda;
    H.turn_dop = d_turn;

    H.dop_raw = dop0;
    H.chi_raw = chi0;

    H.dop_fwd = dop_fwd;
    H.psi_fwd = psi_f;
    H.a1_fwd  = a1_f;
    H.chi_fwd = chi_f;

    H.dop_ret = dop_bwd;
    H.psi_ret = psi_b;
    H.a1_ret  = a1_b;
    H.chi_ret = chi_b;

    H.scan = scan;
    H.opt  = opt;
    H.png  = out_png;
end

% =====================================================================
% pick minimum closest to previous solution (continuation)
% =====================================================================
function pick = pick_closest_min_(mins, prev)
    d = sqrt(sum((mins - prev).^2, 2));
    [~,ix] = min(d);
    pick = mins(ix,:);
end

% =====================================================================
% Loader: recursive, ABS iq/jq, CURRENT chi header (your new format)
%   chi_used = Re(chi) (NO sign flip)
% Header examples:
%   # T = 0
%   # doping = 0.85436393
%   # iq_range = [-133, -132]
%   # jq_range = [-133, -132]
% Data columns:
%   idx iq jq qx qy Re Im [extra...]
% =====================================================================
function D = load_chi_grid_newpatch_absq_recursive_(root_dir, iq_pick, jq_pick)
    root_dir = string(root_dir);
    L = dir(fullfile(root_dir, "**", "chi*.txt"));
    if isempty(L)
        error("No chi*.txt found under: %s", root_dir);
    end

    iq_pick = round(iq_pick);
    jq_pick = round(jq_pick);

    Tvals   = nan(numel(L),1);
    dopvals = nan(numel(L),1);
    chiUsed = nan(numel(L),1);

    n_ok=0; n_fail_header=0; n_fail_read=0; n_fail_match=0; n_oob=0;

    for k = 1:numel(L)
        fpath = fullfile(L(k).folder, L(k).name);

        H = parse_header_newpatch_(fpath);
        if ~H.ok
            n_fail_header = n_fail_header + 1;
            continue;
        end

        if isfinite(H.iq_min) && isfinite(H.iq_max) && isfinite(H.jq_min) && isfinite(H.jq_max)
            if iq_pick < H.iq_min || iq_pick > H.iq_max || jq_pick < H.jq_min || jq_pick > H.jq_max
                n_oob = n_oob + 1;
                continue;
            end
        end

        try
            M = readmatrix(fpath, "FileType","text", "CommentStyle","#");
        catch
            n_fail_read = n_fail_read + 1;
            continue;
        end
        if isempty(M) || size(M,2) < 6
            n_fail_read = n_fail_read + 1;
            continue;
        end

        Cc = detect_cols_(M);

        id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
        if isempty(id)
            n_fail_match = n_fail_match + 1;
            continue;
        end

        ReChi = M(id, Cc.Re);

        Tvals(k)   = H.T_K;
        dopvals(k) = H.dop;
        chiUsed(k) = ReChi;
        n_ok = n_ok + 1;
    end

    mask = isfinite(Tvals) & isfinite(dopvals) & isfinite(chiUsed);
    Tvals   = round(Tvals(mask), 12);
    dopvals = round(dopvals(mask), 12);
    chiUsed = chiUsed(mask);

    fprintf("[load abs] files=%d ok=%d header_fail=%d read_fail=%d oob=%d match_fail=%d valid=%d | iq=%d jq=%d\n", ...
        numel(L), n_ok, n_fail_header, n_fail_read, n_oob, n_fail_match, numel(Tvals), iq_pick, jq_pick);

    if isempty(Tvals)
        D = struct('T_list',zeros(0,1),'dop_list',zeros(0,1),'chi_used_map',[]);
        return;
    end

    T_list   = sort(unique(Tvals));
    dop_list = sort(unique(dopvals));
    NT = numel(T_list);
    Nd = numel(dop_list);

    chi_used_map = nan(NT,Nd);

    for i = 1:numel(Tvals)
        [tfT,iT] = ismember(Tvals(i), T_list);
        [tfD,iD] = ismember(dopvals(i), dop_list);
        if tfT && tfD
            chi_used_map(iT,iD) = chiUsed(i);
        end
    end

    D = struct('T_list',T_list,'dop_list',dop_list,'chi_used_map',chi_used_map);
end

% =====================================================================
% Header parser for your new chi*.txt
% Accepts:
%   # T = ...    or # T_K = ...
%   # doping = ... or # dop = ...
%   # iq_range = [a, b]
%   # jq_range = [a, b]
% =====================================================================
function H = parse_header_newpatch_(fpath)
    H = struct('ok',false,'T_K',NaN,'dop',NaN,'iq_min',NaN,'iq_max',NaN,'jq_min',NaN,'jq_max',NaN);

    fid = fopen(fpath,'r');
    if fid < 0, return; end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";
    inum = "([+\-]?\d+)";

    for t = 1:4000
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        if isempty(ln), continue; end

        ln = strtrim(string(ln));
        if ln == "", continue; end
        if ~startsWith(ln, "#"), break; end

        if ~isfinite(H.T_K)
            tok = regexp(ln, "^\#\s*T\s*=\s*" + num, "tokens","once");
            if isempty(tok)
                tok = regexp(ln, "^\#\s*T_K\s*=\s*" + num, "tokens","once");
            end
            if ~isempty(tok), H.T_K = str2double(tok{1}); end
        end

        if ~isfinite(H.dop)
            tok = regexp(ln, "^\#\s*(doping|dop)\s*=\s*" + num, "tokens","once");
            if ~isempty(tok), H.dop = str2double(tok{2}); end
        end

        if ~isfinite(H.iq_min) || ~isfinite(H.iq_max)
            tok = regexp(ln, "^\#\s*iq_range\s*=\s*\[\s*" + inum + "\s*,\s*" + inum + "\s*\]", "tokens","once");
            if ~isempty(tok)
                H.iq_min = str2double(tok{1});
                H.iq_max = str2double(tok{2});
            end
        end

        if ~isfinite(H.jq_min) || ~isfinite(H.jq_max)
            tok = regexp(ln, "^\#\s*jq_range\s*=\s*\[\s*" + inum + "\s*,\s*" + inum + "\s*\]", "tokens","once");
            if ~isempty(tok)
                H.jq_min = str2double(tok{1});
                H.jq_max = str2double(tok{2});
            end
        end
    end

    H.ok = isfinite(H.T_K) && isfinite(H.dop);
end

% =====================================================================
% Robust column detection
% Supports:
%   idx iq jq qx qy Re Im ...  => iq=2 jq=3 Re=6
%   iq  jq qx qy Re Im ...     => iq=1 jq=2 Re=5
% =====================================================================
function C = detect_cols_(M)
    ncol = size(M,2);
    C = struct('iq',2,'jq',3,'Re',6);

    if ncol < 6
        C.iq = 1; C.jq = 2; C.Re = max(1, min(5,ncol));
        return;
    end

    intlike = @(x) all(abs(x - round(x)) < 1e-9);

    if ncol >= 7
        col1 = M(:,1);
        col2 = M(:,2);
        col3 = M(:,3);

        if intlike(col1) && (min(col1)==0 || min(col1)==1) && intlike(col2) && intlike(col3)
            if any(col2 < 0) || any(col3 < 0) || max(abs(col2)) > 5 || max(abs(col3)) > 5
                C.iq = 2; C.jq = 3; C.Re = 6;
                return;
            end
        end
    end

    C.iq = 1;
    C.jq = 2;
    C.Re = 5;
end

% =====================================================================
% dense scan -> seeds -> fminsearch -> cluster -> Hessian PD
% =====================================================================
function mins = find_local_minima_2d_denseSeeds_(F, psi1_lim, psi2_lim, scan, opt)
    psi1 = linspace(psi1_lim(1), psi1_lim(2), scan.Npsi1);
    psi2 = linspace(psi2_lim(1), psi2_lim(2), scan.Npsi2);
    [P1,P2] = meshgrid(psi1, psi2);

    Fgrid = arrayfun(@(u,v) safe_F_(F,u,v,psi1_lim,psi2_lim), P1, P2);

    if scan.use_gaussian_smooth
        Fgrid = gaussian_smooth_2d_(Fgrid, scan.smooth_sigma);
    end

    mask_min = discrete_local_min_mask_(Fgrid);
    [rr,cc] = find(mask_min);

    if isempty(rr)
        seeds = [0,0];
    else
        seeds = [P1(sub2ind(size(P1), rr, cc)), P2(sub2ind(size(P2), rr, cc))];
        E = Fgrid(sub2ind(size(Fgrid), rr, cc));
        [~,ix] = sort(E,'ascend');
        ix = ix(1:min(scan.keep_topK, numel(ix)));
        seeds = seeds(ix,:);
    end

    if opt.force_origin_and_nearby
        d = opt.nearby_delta;
        seeds = [seeds;
                 0,0;
                 d,0; -d,0; 0,d; 0,-d];
    end

    obj = @(x) safe_F_(F, x(1), x(2), psi1_lim, psi2_lim);

    fopt = optimset('Display','off', ...
        'MaxIter', opt.max_iter, ...
        'TolFun', opt.tol_fun, ...
        'TolX',   opt.tol_x);

    cand = zeros(0,2);
    for k = 1:size(seeds,1)
        x0 = seeds(k,:);
        for rep = 1:opt.seed_repeats
            xstart = x0 + opt.seed_jitter * randn(1,2);
            xstart(1) = min(max(xstart(1), psi1_lim(1)), psi1_lim(2));
            xstart(2) = min(max(xstart(2), psi2_lim(1)), psi2_lim(2));
            try
                xhat = fminsearch(obj, xstart, fopt);
            catch
                continue;
            end
            xhat(1) = min(max(xhat(1), psi1_lim(1)), psi1_lim(2));
            xhat(2) = min(max(xhat(2), psi2_lim(1)), psi2_lim(2));
            cand(end+1,:) = xhat; %#ok<AGROW>
        end
    end

    if isempty(cand)
        mins = zeros(0,2);
        return;
    end

    uniq = cluster_points_(cand, opt.cluster_tol);

    keep = false(size(uniq,1),1);
    for i = 1:size(uniq,1)
        x = uniq(i,:);
        H = hessian_fd_(@(u,v) F(u,v), x(1), x(2), opt.fd_h);
        ev = eig((H+H')/2);
        if all(real(ev) > opt.min_eig_eps)
            keep(i) = true;
        end
    end

    mins = uniq(keep,:);
end

function val = safe_F_(F, psi1, psi2, lim1, lim2)
    if psi1 < lim1(1) || psi1 > lim1(2) || psi2 < lim2(1) || psi2 > lim2(2)
        val = 1e30 + 1e10*(abs(psi1)+abs(psi2));
        return;
    end
    val = F(psi1, psi2);
    if ~isfinite(val), val = 1e30; end
end

function uniq = cluster_points_(P, tol)
    uniq = zeros(0,2);
    for i = 1:size(P,1)
        x = P(i,:);
        if isempty(uniq)
            uniq = x;
        else
            d = sqrt(sum((uniq - x).^2, 2));
            if all(d > tol)
                uniq(end+1,:) = x; %#ok<AGROW>
            end
        end
    end
end

function H = hessian_fd_(f, x, y, h)
    f00 = f(x,y);
    fxx = (f(x+h,y) - 2*f00 + f(x-h,y))/h^2;
    fyy = (f(x,y+h) - 2*f00 + f(x,y-h))/h^2;
    fxy = (f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)) / (4*h^2);
    H = [fxx, fxy; fxy, fyy];
end

function mask = discrete_local_min_mask_(A)
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

function B = gaussian_smooth_2d_(A, sigma)
    if sigma <= 0, B = A; return; end
    rad = max(1, ceil(3*sigma));
    x = (-rad:rad);
    g = exp(-(x.^2)/(2*sigma^2));
    g = g / sum(g);
    B = conv2(A, g, 'same');
    B = conv2(B, g', 'same');
end
