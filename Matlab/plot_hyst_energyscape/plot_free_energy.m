function out = plot_free_energy()
% plot_free_energy
%
% Read:
%   T_target, lambda, iq/jq, psi limits from make_realchi_params(true)
%
% Manual input:
%   doping_target
%
% chi_used is read from chi*.txt recursively by finding the nearest
% (T, doping) point to (par.T_target, doping_target).
%
% Plot 3 figures:
%   1) total free-energy heatmap F(psi1,psi2)
%   2) electronic free energy F_psi(psi1)
%   3) lattice free energy F_X(psi2)

    % =========================================================
    % 1. params + manual input
    % =========================================================
    par = make_realchi_params(true);

    T_target = par.T_target;
    lambda   = par.lambda;

    % ---- only this is manually chosen ----
    doping_target = -0.15;

    iq_pick = round(par.iq_pick);
    jq_pick = round(par.jq_pick);

    % nearest-point metric weights
    T_scale   = 0.5;
    dop_scale = 0.05;

    max_abs_dT   = inf;
    max_abs_ddop = inf;

    % plotting resolution
    N1 = 201;
    N2 = 201;
    fs = par.plot.fontSize;

    find_minimum = true;

    % =========================================================
    % 2. choose data folder
    % =========================================================
    default_root = par.io.default_root;
    if ~isfolder(default_root)
        default_root = string(pwd);
    end

    root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
    if isequal(root,0)
        error('User cancelled.');
    end
    root = string(root);

    % =========================================================
    % 3. load scattered chi(T,doping)
    % =========================================================
    S = load_chi_scattered_T_doping_(root, iq_pick, jq_pick);

    if isempty(S.T)
        error('No valid chi points found.');
    end

    % =========================================================
    % 4. find nearest data point
    % =========================================================
    dT   = S.T - T_target;
    ddop = S.doping - doping_target;

    metric = (dT ./ T_scale).^2 + (ddop ./ dop_scale).^2;
    [metric_min, idx_best] = min(metric);

    T_used      = S.T(idx_best);
    doping_used = S.doping(idx_best);
    chi_used    = S.chi(idx_best);
    file_used   = S.file(idx_best);

    abs_dT   = abs(T_used - T_target);
    abs_ddop = abs(doping_used - doping_target);

    if abs_dT > max_abs_dT || abs_ddop > max_abs_ddop
        error('Nearest point is too far: dT=%.6g, ddop=%.6g', abs_dT, abs_ddop);
    end

    fprintf('\n========== nearest chi point ==========\n');
    fprintf('Requested T       = %.12g\n', T_target);
    fprintf('Requested doping  = %.12g\n', doping_target);
    fprintf('Matched   T       = %.12g\n', T_used);
    fprintf('Matched   doping  = %.12g\n', doping_used);
    fprintf('chi_used          = %.12g\n', chi_used);
    fprintf('lambda            = %.12g\n', lambda);
    fprintf('metric            = %.12g\n', metric_min);
    fprintf('source file       = %s\n', file_used);
    fprintf('=======================================\n\n');

    % =========================================================
    % 5. build coeffs
    % =========================================================
    coef = make_realchi_coeff(par);
    C = coef.eval(T_used, doping_used, chi_used);

    a1 = C.a1;
    b1 = C.b1;
    a2 = C.a2;
    b2 = C.b2;
    c2 = C.c2;

    delta = (5/6) * b2^2 - a2 * c2;

    % =========================================================
    % 6. define free energies
    % =========================================================
    Fpsi = @(psi1) ...
        0.5 * a1 * psi1.^2 + (1/factorial(4)) * b1 * psi1.^4;

    FX = @(psi2) ...
        0.5 * a2 * psi2.^2 ...
        + (1/factorial(4)) * b2 * psi2.^4 ...
        + (1/factorial(6)) * c2 * psi2.^6;

    F = @(psi1, psi2) ...
        0.5 * a1 * psi1.^2 ...
        + (1/factorial(4)) * b1 * psi1.^4 ...
        + 0.5 * a2 * psi2.^2 ...
        + (1/factorial(4)) * b2 * psi2.^4 ...
        + (1/factorial(6)) * c2 * psi2.^6 ...
        + lambda * psi1 .* psi2;

    grad_F = @(psi1, psi2) [ ...
        a1 * psi1 + (1/6) * b1 * psi1.^3 + lambda * psi2; ...
        a2 * psi2 + (1/6) * b2 * psi2.^3 + (1/120) * c2 * psi2.^5 + lambda * psi1];

    % =========================================================
    % 7. grids
    % =========================================================
    psi1_vals = linspace(par.psi1_lim(1), par.psi1_lim(2), N1);
    psi2_vals = linspace(par.psi2_lim(1), par.psi2_lim(2), N2);
    [P1, P2] = meshgrid(psi1_vals, psi2_vals);
    Fv = arrayfun(F, P1, P2);

    % =========================================================
    % 8. optional minimum search
    % =========================================================
    psi1_star = NaN;
    psi2_star = NaN;
    F_star    = NaN;

    if find_minimum
        obj = @(x) F(x(1), x(2));

        x0_list = [
             0,  0;
             1,  1;
            -1,  1;
             1, -1;
            -1, -1;
             5,  5;
            -5,  5;
             5, -5;
            -5, -5
        ];

        best_x = [0, 0];
        best_f = inf;
        opts = optimset('Display', 'off');

        for k = 1:size(x0_list,1)
            try
                xk = fminsearch(obj, x0_list(k,:), opts);
                fk = obj(xk);
                if isfinite(fk) && fk < best_f
                    best_f = fk;
                    best_x = xk;
                end
            catch
            end
        end

        psi1_star = best_x(1);
        psi2_star = best_x(2);
        F_star    = best_f;
    end

    % =========================================================
    % 9. figure 1: total free-energy heatmap
    % =========================================================
    fig1 = figure('Color','w','Name','Free Energy Heatmap');
    ax1 = axes(fig1);

    pcolor(ax1, P1, P2, Fv);
    shading(ax1, 'interp');
    colormap(ax1, 'turbo');
    colorbar(ax1, 'Location', 'eastoutside');

    hold(ax1, 'on');
    contour(ax1, P1, P2, Fv, 24, 'k-', 'LineWidth', 0.5);
    if isfinite(psi1_star) && isfinite(psi2_star)
        plot(ax1, psi1_star, psi2_star, 'ro', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 7);
    end
    hold(ax1, 'off');

    xlim(ax1, par.psi1_lim);
    ylim(ax1, par.psi2_lim);
    xlabel(ax1, '\psi_1');
    ylabel(ax1, '\psi_2');
    set(ax1, 'FontSize', fs, 'LineWidth', 1.2, 'Box', 'on');

    title(ax1, sprintf('F(\\psi_1,\\psi_2): T=%.3f, doping=%.3f, \\lambda=%.3f', ...
        T_used, doping_used, lambda), 'FontWeight', 'normal');

    % =========================================================
    % 10. figure 2: electronic free energy
    % =========================================================
    fig2 = figure('Color','w','Name','Electronic Free Energy');
    ax2 = axes(fig2);

    y1 = Fpsi(psi1_vals);
    plot(ax2, psi1_vals, y1, 'LineWidth', 2);
    grid(ax2, 'on');
    box(ax2, 'on');
    xlim(ax2, par.psi1_lim);
    xlabel(ax2, '\psi_1');
    ylabel(ax2, 'F_\psi(\psi_1)');
    set(ax2, 'FontSize', fs, 'LineWidth', 1.2);

    [~, idx1] = min(y1);
    hold(ax2, 'on');
    plot(ax2, psi1_vals(idx1), y1(idx1), 'ro', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 6);
    hold(ax2, 'off');

    title(ax2, sprintf('Electronic: a1=%.6g, b1=%.6g, chi=%.6g', a1, b1, chi_used), ...
        'FontWeight', 'normal');

    % =========================================================
    % 11. figure 3: lattice free energy
    % =========================================================
    fig3 = figure('Color','w','Name','Lattice Free Energy');
    ax3 = axes(fig3);

    y2 = FX(psi2_vals);
    plot(ax3, psi2_vals, y2, 'LineWidth', 2);
    grid(ax3, 'on');
    box(ax3, 'on');
    xlim(ax3, par.psi2_lim);
    xlabel(ax3, '\psi_2');
    ylabel(ax3, 'F_X(\psi_2)');
    set(ax3, 'FontSize', fs, 'LineWidth', 1.2);

    [~, idx2] = min(y2);
    hold(ax3, 'on');
    plot(ax3, psi2_vals(idx2), y2(idx2), 'ro', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 6);
    hold(ax3, 'off');

    title(ax3, sprintf('Lattice: a2=%.6g, b2=%.6g, c2=%.6g, \\Delta=%.6g', ...
        a2, b2, c2, delta), 'FontWeight', 'normal');

    % =========================================================
    % 12. output
    % =========================================================
    out = struct();
    out.T_target = T_target;
    out.doping_target = doping_target;
    out.T_used = T_used;
    out.doping_used = doping_used;
    out.lambda = lambda;
    out.chi_used = chi_used;
    out.file_used = file_used;
    out.a1 = a1;
    out.b1 = b1;
    out.a2 = a2;
    out.b2 = b2;
    out.c2 = c2;
    out.delta = delta;
    out.F = F;
    out.Fpsi = Fpsi;
    out.FX = FX;
    out.grad_F = grad_F;
    out.psi1_star = psi1_star;
    out.psi2_star = psi2_star;
    out.F_star = F_star;
end

% =====================================================================
% helper: load scattered chi(T,doping)
% =====================================================================
function S = load_chi_scattered_T_doping_(root_dir, iq_pick, jq_pick)

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

    n_ok = 0;

    for k = 1:n
        fpath = string(fullfile(L(k).folder, L(k).name));

        H = parse_header_T_and_doping_robust_(fpath);
        if ~H.ok
            continue;
        end

        try
            M = read_chi_numeric_skiphash_(fpath);
        catch
            continue;
        end

        if isempty(M) || size(M,2) < 5
            continue;
        end

        Cc = detect_cols_(M);
        id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
        if isempty(id)
            continue;
        end

        Tvals(k) = H.T_K;
        doping_header(k) = H.doping;
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

    fprintf('[load chi] files=%d, valid=%d\n', n, n_ok);
end

% =====================================================================
% helper: parse header
% =====================================================================
function H = parse_header_T_and_doping_robust_(fpath)

    H = struct('ok',false,'T_K',NaN,'doping',NaN);

    fid = fopen(fpath,'r');
    if fid < 0
        return;
    end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';

    while true
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
            tok = regexp(s, ['T\s*=\s*' num], 'tokens','once');
            if isempty(tok)
                tok = regexp(s, ['T_K\s*=\s*' num], 'tokens','once');
            end
            if ~isempty(tok)
                H.T_K = str2double(tok{1});
            end
        end

        if ~isfinite(H.doping)
            tok = regexp(s, ['doping\s*=\s*' num], 'tokens','once');
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

% =====================================================================
% helper: read numeric body
% =====================================================================
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
end

% =====================================================================
% helper: detect columns
% =====================================================================
function C = detect_cols_(M)

    ncol = size(M,2);

    if ncol >= 6
        C.iq = 2;
        C.jq = 3;
        C.Re = 6;
    else
        C.iq = 1;
        C.jq = 2;
        C.Re = min(5, ncol);
    end
end