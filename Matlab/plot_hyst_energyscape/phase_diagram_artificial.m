function out = phase_diagram_artificial()
% plot_phase_diagram_free_energy
% Phase diagram in (n, T) plane by counting local minima of full 2D free energy
% F(psi, X) = free_energy(a1,b1,a2,b2,c2,lambda).
%
% 关键改动：
%   1) 对每个 (T,n) 先在更密的 (psi,X) 网格上离散采样 F；
%   2) 从离散网格里抓“局部极小值格点”作为 seeds；
%   3) 对这些 seeds 做 fminsearch (带 jitter 多次)；
%   4) cluster 去重 + Hessian 正定筛选，得到最终极小值集合。
%
% Output:
%   out.phase_map (NT x Nn) integer in {0,1,2,3,4}
%   out.nmin_map  (NT x Nn) #minima for debug
%   out.T_list, out.n_list

    % =====================================================================
    % (0) ALL TUNABLE PARAMETERS (put everything here)
    % =====================================================================

    % -------- phase diagram grid (T,n) --------
    par = struct();
    par.T_list = linspace(50, 90, 41);     % T sampling for phase diagram
    par.n_list = linspace(-10, 10, 81);     % n sampling for phase diagram
    par.lambda = 0.5;                       % fix lambda for this phase diagram

    % -------- search domain in (psi,X) --------
    par.psi1_lim = [-12, 12];               % psi range to search minima
    par.psi2_lim = [-8,  8];                % X   range to search minima

    % -------- dense discrete grid for scanning F(psi,X) --------
    par.scan_Npsi1 = 81;                    % dense grid points along psi
    par.scan_Npsi2 = 81;                    % dense grid points along X
    par.scan_use_gaussian_smooth = true;    % 先对网格F做轻微平滑再找局部最小(抗噪/抗离散抖动)
    par.scan_smooth_sigma = 0.7;            % 平滑强度(单位：grid step)，0.5~1.0 常用
    par.scan_keep_topK = 30;                % 只保留能量最低的前K个离散局部极小点当 seeds（控成本）

    % -------- minimization settings (fminsearch) --------
    opt = struct();
    opt.seed_jitter  = 0.12;                % seed around each discrete min (Gaussian jitter)
    opt.seed_repeats = 2;                   % repeats per seed (with different jitter)
    opt.max_iter     = 400;                 % fminsearch max iterations
    opt.tol_fun      = 1e-10;               % fminsearch TolFun
    opt.tol_x        = 1e-8;                % fminsearch TolX

    % -------- dedup / acceptance --------
    opt.cluster_tol  = 3e-2;                % dedup tolerance in (psi,X)
    opt.fd_h         = 2e-3;                % finite diff step for Hessian check
    opt.min_eig_eps  = 1e-6;                % Hessian eigenvalues must be > this
    opt.origin_tol   = 5e-2;                % classify "origin minimum" if within this distance

    % -------- extra seeds always included --------
    opt.force_origin_and_nearby = true;     % always include origin and small neighborhood
    opt.nearby_delta = 0.2;                 % nearby seed offset

    % -------- coefficient formulas --------
    beta = @(T) 2 - (T-70)/20;
    n_1  = @(n) (n - abs(n))/2 + (abs(n) + n)/2; % == n

    a1_formula = @(T, n) -2.5 * beta(T) - 1.7*(n_1(n-2)+1.3) + 2;
    b1_formula = @(T) 1;
    b2_formula = @(T) -1 * beta(T);
    a2_formula = @(T,n) 0.5 * (beta(T).^2) + 1.4 - 0.075*(n_1(n-2)+1.3);
    c2_formula = @(T) 1;

    % -------- plotting --------
    par.plot_colormap = turbo;              % colormap for phase map
    par.overlay_boundary = true;            % overlay boundary dots

    % =====================================================================
    % (A) Allocate results
    % =====================================================================
    T_list = par.T_list;
    n_list = par.n_list;

    NT = numel(T_list);
    Nn = numel(n_list);

    phase_map = zeros(NT, Nn); % 0 means unknown/unclassified
    nmin_map  = zeros(NT, Nn);

    fprintf('Scanning phase diagram: NT=%d, Nn=%d (lambda=%.3g)\n', NT, Nn, par.lambda);

    % =====================================================================
    % (B) Main loop
    % =====================================================================
    for iT = 1:NT
        T = T_list(iT);

        for in = 1:Nn
            n = n_list(in);

            a1 = a1_formula(T, n);
            b1 = b1_formula(T);
            a2 = a2_formula(T, n);
            b2 = b2_formula(T);
            c2 = c2_formula(T);

            F = free_energy(a1, b1, a2, b2, c2, par.lambda); % F(psi,X)

            mins = find_local_minima_2d_denseSeeds(F, par.psi1_lim, par.psi2_lim, par, opt);

            nmin_map(iT,in)  = size(mins,1);
            phase_map(iT,in) = classify_by_minima(mins, opt);
        end

        fprintf('  done T=%g (%d/%d)\n', T, iT, NT);
    end

    % =====================================================================
    % (C) Plot phase diagram
    % =====================================================================
    fig = figure('Name','Phase diagram (n vs T)','Color','w'); %#ok<NASGU>
    imagesc(n_list, T_list, phase_map); axis xy;
    xlabel('n'); ylabel('T');
    title(sprintf('Phase diagram by counting local minima (\\lambda=%.3g)', par.lambda));
    colormap(par.plot_colormap);
    cb = colorbar;
    cb.Ticks = 0:4;
    cb.TickLabels = {'0:unknown','1','2','3','4'};

    if par.overlay_boundary
        hold on;
        B = boundary_mask(phase_map);
        [r,c] = find(B);
        plot(n_list(c), T_list(r), 'k.', 'MarkerSize', 4);
        hold off;
    end

    out = struct();
    out.phase_map = phase_map;
    out.nmin_map  = nmin_map;
    out.T_list    = T_list;
    out.n_list    = n_list;
    out.lambda    = par.lambda;
end

% =====================================================================
% helper: dense grid scan -> pick discrete local minima as seeds -> fminsearch
% =====================================================================
function mins = find_local_minima_2d_denseSeeds(F, psi1_lim, psi2_lim, par, opt)

    % ---------- (1) dense discrete scan ----------
    psi1 = linspace(psi1_lim(1), psi1_lim(2), par.scan_Npsi1);
    psi2 = linspace(psi2_lim(1), psi2_lim(2), par.scan_Npsi2);
    [P1,P2] = meshgrid(psi1, psi2);

    % evaluate F on grid (vectorized with arrayfun; keeps compatibility)
    Fgrid = arrayfun(@(u,v) safe_F(F,u,v,psi1_lim,psi2_lim), P1, P2);

    % optional smoothing (reduce spurious "minima" caused by discretization)
    if par.scan_use_gaussian_smooth
        Fgrid = gaussian_smooth_2d(Fgrid, par.scan_smooth_sigma);
    end

    % ---------- (2) find discrete local minima on grid ----------
    mask_min = discrete_local_min_mask(Fgrid);

    [rr,cc] = find(mask_min);
    if isempty(rr)
        % fallback: at least try origin
        seeds = [0,0];
    else
        seeds = [P1(sub2ind(size(P1), rr, cc)), P2(sub2ind(size(P2), rr, cc))];

        % keep topK lowest-energy discrete minima (control runtime)
        E = Fgrid(sub2ind(size(Fgrid), rr, cc));
        [~,ix] = sort(E, 'ascend');
        ix = ix(1:min(par.scan_keep_topK, numel(ix)));
        seeds = seeds(ix,:);
    end

    % always include origin & nearby points (helps when scan misses origin basin)
    if opt.force_origin_and_nearby
        d = opt.nearby_delta;
        seeds = [seeds;
                 0,0;
                 d,0; -d,0; 0,d; 0,-d];
    end

    % ---------- (3) run fminsearch from seeds ----------
    obj = @(x) safe_F(F, x(1), x(2), psi1_lim, psi2_lim);

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

    % ---------- (4) cluster (dedup) ----------
    uniq = cluster_points(cand, opt.cluster_tol);

    % ---------- (5) keep only true minima by Hessian positive definite ----------
    keep = false(size(uniq,1),1);
    for i = 1:size(uniq,1)
        x = uniq(i,:);
        H = hessian_fd(@(u,v) F(u,v), x(1), x(2), opt.fd_h);
        ev = eig((H+H')/2);
        if all(real(ev) > opt.min_eig_eps)
            keep(i) = true;
        end
    end

    mins = uniq(keep,:);
end

% =====================================================================
% utilities
% =====================================================================
function val = safe_F(F, psi1, psi2, lim1, lim2)
    if psi1 < lim1(1) || psi1 > lim1(2) || psi2 < lim2(1) || psi2 > lim2(2)
        val = 1e30 + 1e10*(abs(psi1)+abs(psi2));
        return;
    end
    val = F(psi1, psi2);
    if ~isfinite(val)
        val = 1e30;
    end
end

function uniq = cluster_points(P, tol)
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

function H = hessian_fd(f, x, y, h)
    f00 = f(x,y);
    fxx = (f(x+h,y) - 2*f00 + f(x-h,y))/h^2;
    fyy = (f(x,y+h) - 2*f00 + f(x,y-h))/h^2;
    fxy = (f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)) / (4*h^2);
    H = [fxx, fxy; fxy, fyy];
end

function mask = discrete_local_min_mask(A)
% Discrete local minima mask on 2D grid:
% A(i,j) is a local minimum if it is <= all its 8 neighbors.
% Edge points are excluded to keep logic simple (can be extended if needed).
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

function B = gaussian_smooth_2d(A, sigma)
% Simple Gaussian smoothing using separable convolution.
% sigma is in units of grid steps.
    if sigma <= 0
        B = A; return;
    end

    rad = max(1, ceil(3*sigma));
    x = (-rad:rad);
    g = exp(-(x.^2)/(2*sigma^2));
    g = g / sum(g);

    % convolve rows then cols
    B = conv2(A, g, 'same');
    B = conv2(B, g', 'same');
end

% =====================================================================
% classify phases by number of minima + origin-min check
% =====================================================================
function phase = classify_by_minima(mins, opt)
    nmin = size(mins,1);

    if nmin == 0
        phase = 0;
        return;
    end

    d0 = sqrt(sum(mins.^2,2));
    has_origin_min = any(d0 < opt.origin_tol);

    switch nmin
        case 1
            if has_origin_min
                phase = 1;
            else
                phase = 0;
            end
        case 2
            phase = 4;
        case 3
            if has_origin_min
                phase = 2;
            else
                phase = 0;
            end
        case 4
            phase = 3;
        otherwise
            phase = 0;
    end
end

function B = boundary_mask(M)
    [R,C] = size(M);
    B = false(R,C);
    for r = 1:R
        for c = 1:C
            v = M(r,c);
            if r>1   && M(r-1,c)~=v, B(r,c)=true; end
            if r<R   && M(r+1,c)~=v, B(r,c)=true; end
            if c>1   && M(r,c-1)~=v, B(r,c)=true; end
            if c<C   && M(r,c+1)~=v, B(r,c)=true; end
        end
    end
end
