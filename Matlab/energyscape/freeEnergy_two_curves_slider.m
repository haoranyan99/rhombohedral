function freeEnergy_two_curves_slider()
% freeEnergy_two_curves_fixedT_alpha_slider
%
% UI for plotting:
%   1) electronic free energy F_psi(psi1)
%   2) lattice free energy    F_X(psi2)
%
% Controls:
%   - fixed T input
%   - alpha slider
%   - doping(folder) or mu(folder) slider
%   - iq/jq selection
%
% Requires:
%   - make_realchi_params
%   - make_realchi_coeff
%
% Data source:
%   choose a folder containing chi*.txt recursively

    % =========================================================
    % 1. base params
    % =========================================================
    par_base = make_realchi_params(true);

    % ---- UI defaults ----
    T_fixed_default = 2;

    alpha_min = 0.01;
    alpha_max = max(5.0, 2.0 * par_base.Lat_alpha);
    alpha0    = par_base.Lat_alpha;

    % ---- default root ----
    default_root = "/Users/haoranyan/rg_master/data/";
    if isfield(par_base,'io') && isfield(par_base.io,'default_root') ...
            && isfolder(par_base.io.default_root)
        default_root = string(par_base.io.default_root);
    end
    if ~isfolder(default_root)
        default_root = string(pwd);
    end

    % =========================================================
    % 2. choose folder
    % =========================================================
    root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
    if isequal(root,0)
        return;
    end
    root = string(root);

    % =========================================================
    % 3. load chi grid
    % =========================================================
    G = load_chi_grid_folderUI_(root, par_base.iq_pick, par_base.jq_pick);

    T_list  = G.T_list;
    U_list  = G.U_list;
    chi_map = G.chi_map;
    dop_map = G.doping_map;
    muf_map = G.mu_folder_map;
    u_tag   = G.u_tag;

    if isempty(T_list) || isempty(U_list)
        error('No valid chi grid found under selected folder.');
    end

    [~, it0] = min(abs(T_list - T_fixed_default));
    iu0 = nearest_valid_in_row_(chi_map, it0, 1);
    if ~isfinite(chi_map(it0, iu0))
        [it0, iu0] = find_first_valid_(chi_map);
    end

    % =========================================================
    % 4. UI
    % =========================================================
    fig = uifigure( ...
        'Name', sprintf('Two Free Energies | fixed T + alpha slider (%s)', u_tag), ...
        'Position', [80 60 1220 720]);

    % ---- axes ----
    axPsi = uiaxes(fig, 'Position', [60 300 480 340], 'FontSize', 15);
    axPsi.TickLabelInterpreter = 'latex';
    xlabel(axPsi, '$\psi_1$', 'Interpreter', 'latex');
    ylabel(axPsi, '$F_{\psi}(\psi_1)$', 'Interpreter', 'latex');
    title(axPsi, '$F_{\psi}=\frac12 a_1\psi_1^2+\frac{1}{4!}b_1\psi_1^4$', ...
        'Interpreter', 'latex', 'FontWeight', 'normal');

    axX = uiaxes(fig, 'Position', [660 300 480 340], 'FontSize', 15);
    axX.TickLabelInterpreter = 'latex';
    xlabel(axX, '$\psi_2$', 'Interpreter', 'latex');
    ylabel(axX, '$F_{X}(\psi_2)$', 'Interpreter', 'latex');
    title(axX, '$F_X=\frac12 a_2\psi_2^2+\frac{1}{4!}b_2\psi_2^4+\frac{1}{6!}c_2\psi_2^6$', ...
        'Interpreter', 'latex', 'FontWeight', 'normal');

    % ---- parameter display ----
    coeff_box = uitextarea(fig, ...
        'Position', [60 25 1080 240], ...
        'Editable', 'off', ...
        'FontSize', 15);

    % =========================================================
    % 5. controls
    % =========================================================
    y_top = 255;
    dy    = 55;

    % ---- fixed T input ----
    uilabel(fig, ...
        'Position', [60 y_top+dy-10 90 24], ...
        'Text', 'T fixed', ...
        'FontSize', 15);

    T_input = uieditfield(fig, 'numeric', ...
        'Position', [130 y_top+dy-12 120 30], ...
        'Value', T_list(it0), ...
        'FontSize', 15);

    T_text = uilabel(fig, ...
        'Position', [270 y_top+dy-12 340 24], ...
        'Text', '', ...
        'FontSize', 15);

    % ---- alpha slider ----
    uilabel(fig, ...
        'Position', [60 y_top-10 90 24], ...
        'Text', '\alpha', ...
        'FontSize', 15);

    alpha_slider = uislider(fig, ...
        'Position', [130 y_top 380 3], ...
        'Limits', [alpha_min, alpha_max], ...
        'Value', alpha0);
    alpha_slider.MajorTicks = [];
    alpha_slider.MinorTicks = [];

    alpha_text = uilabel(fig, ...
        'Position', [530 y_top-14 280 24], ...
        'Text', '', ...
        'FontSize', 15);

    % ---- U slider ----
    uilabel(fig, ...
        'Position', [60 y_top-dy-10 90 24], ...
        'Text', char(u_tag), ...
        'FontSize', 15);

    U_slider = uislider(fig, ...
        'Position', [130 y_top-dy 380 3], ...
        'Limits', [1, max(1, numel(U_list))], ...
        'Value', iu0);
    U_slider.MajorTicks = [];
    U_slider.MinorTicks = [];

    U_text = uilabel(fig, ...
        'Position', [530 y_top-dy-14 360 24], ...
        'Text', '', ...
        'FontSize', 15);

    % ---- iq jq ----
    uilabel(fig, 'Position', [930 255 30 22], 'Text', 'iq', 'FontSize', 15);
    iq_input = uieditfield(fig, 'numeric', ...
        'Position', [920 225 80 30], ...
        'Value', par_base.iq_pick, ...
        'FontSize', 15);

    uilabel(fig, 'Position', [1020 255 30 22], 'Text', 'jq', 'FontSize', 15);
    jq_input = uieditfield(fig, 'numeric', ...
        'Position', [1010 225 80 30], ...
        'Value', par_base.jq_pick, ...
        'FontSize', 15);

    apply_button = uibutton(fig, 'push', ...
        'Position', [920 180 170 34], ...
        'Text', 'Apply iq/jq', ...
        'FontSize', 15);

    % =========================================================
    % 6. callbacks
    % =========================================================
    T_input.ValueChangedFcn       = @(~,~) update_();
    alpha_slider.ValueChangedFcn  = @(~,~) update_();
    U_slider.ValueChangedFcn      = @(~,~) update_();
    apply_button.ButtonPushedFcn  = @(~,~) apply_iqjq_();

    refresh_labels_();
    update_();

    % =========================================================
    % nested functions
    % =========================================================
    function refresh_labels_()
        iu = clamp_(round(U_slider.Value), numel(U_list));
        U_slider.Value = iu;

        [~, it_show] = min(abs(T_list - T_input.Value));
        T_text.Text = sprintf('nearest available T = %.6g K', T_list(it_show));

        alpha_text.Text = sprintf('\\alpha = %.6g', alpha_slider.Value);

        uval = U_list(iu);
        if strcmpi(u_tag, 'doping')
            U_text.Text = sprintf('doping(folder) = %.6g', uval);
        else
            U_text.Text = sprintf('mu(folder) = %.6g', uval);
        end
    end

    function update_()
        refresh_labels_();

        [~, it] = min(abs(T_list - T_input.Value));
        iu = clamp_(round(U_slider.Value), numel(U_list));

        if ~isfinite(chi_map(it, iu))
            iu = nearest_valid_in_row_(chi_map, it, iu);
            U_slider.Value = iu;
            refresh_labels_();

            if ~isfinite(chi_map(it, iu))
                render_empty_(sprintf('No valid point at T = %.6g', T_list(it)));
                return;
            end
        end

        T_now = T_list(it);
        u_now = U_list(iu);
        chi_used = chi_map(it, iu);
        dop_eval = dop_map(it, iu);
        mu_folder_value = muf_map(it, iu);
        alpha_now = alpha_slider.Value;

        par_local = par_base;
        par_local.Lat_alpha = alpha_now;
        coef_local = make_realchi_coeff(par_local);

        C = coef_local.eval(T_now, dop_eval, chi_used);

        draw_curves_(par_local, T_now, u_now, chi_used, dop_eval, mu_folder_value, C, alpha_now);
    end

    function draw_curves_(par_local, T_now, u_now, chi_used, dop_eval, mu_folder_value, C, alpha_now)
        cla(axPsi);
        cla(axX);

        % ---- free energies ----
        Fpsi = @(psi) 0.5 * C.a1 * psi.^2 + (1/factorial(4)) * C.b1 * psi.^4;
        FX   = @(X)   0.5 * C.a2 * X.^2 ...
                    + (1/factorial(4)) * C.b2 * X.^4 ...
                    + (1/factorial(6)) * C.c2 * X.^6;

        % ---- x ranges ----
        if isfield(par_local, 'psi1_lim') && numel(par_local.psi1_lim)==2
            x1 = linspace(par_local.psi1_lim(1), par_local.psi1_lim(2), 2001);
            xlim(axPsi, par_local.psi1_lim);
        else
            x1 = linspace(-10, 10, 2001);
            xlim(axPsi, [x1(1), x1(end)]);
        end

        if isfield(par_local, 'psi2_lim') && numel(par_local.psi2_lim)==2
            x2 = linspace(par_local.psi2_lim(1), par_local.psi2_lim(2), 4001);
            xlim(axX, par_local.psi2_lim);
        else
            x2 = linspace(-20, 20, 4001);
            xlim(axX, [x2(1), x2(end)]);
        end

        y1 = Fpsi(x1);
        y2 = FX(x2);

        % ---- electronic plot ----
        plot(axPsi, x1, y1, 'LineWidth', 2);
        grid(axPsi, 'on');
        box(axPsi, 'on');
        [~, idx1] = min(y1);
        hold(axPsi, 'on');
        plot(axPsi, x1(idx1), y1(idx1), 'o', ...
            'MarkerSize', 7, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
        text(axPsi, x1(idx1), y1(idx1), ...
            sprintf('  min: (%.4g, %.4g)', x1(idx1), y1(idx1)), ...
            'FontSize', 12, 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
        hold(axPsi, 'off');

        % ---- lattice plot with ALL local minima ----
        plot(axX, x2, y2, 'LineWidth', 2);
        grid(axX, 'on');
        box(axX, 'on');
        hold(axX, 'on');

        imin = find_local_minima_(y2);
        if isempty(imin)
            [~, idx0] = min(y2);
            imin = idx0;
        end

        min_pos = x2(imin);
        min_val = y2(imin);

        % sort by position
        [min_pos, order] = sort(min_pos);
        min_val = min_val(order);
        imin = imin(order);

        plot(axX, min_pos, min_val, 'o', ...
            'MarkerSize', 7, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

        for k = 1:numel(min_pos)
            text(axX, min_pos(k), min_val(k), ...
                sprintf('  min%d: (%.4g, %.4g)', k, min_pos(k), min_val(k)), ...
                'FontSize', 12, 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
        end

        hold(axX, 'off');

        % ---- parameter box ----
        delta_now = (5/6) * C.b2^2 - C.a2 * C.c2;

        lines = {
            sprintf('iq = %d,   jq = %d', par_local.iq_pick, par_local.jq_pick)
            sprintf('T = %.12g   (requested %.12g)', T_now, T_input.Value)
            sprintf('alpha = %.12g', alpha_now)
        };

        if strcmpi(u_tag, 'doping')
            lines{end+1} = sprintf('doping(folder) = %.12g', u_now);
        else
            lines{end+1} = sprintf('mu(folder) = %.12g', u_now);
        end

        lines = [lines; {
            sprintf('doping(header) = %.12g', dop_eval)
            sprintf('chi = %.12g', chi_used)
            sprintf('a1 = %.12g', C.a1)
            sprintf('b1 = %.12g', C.b1)
            sprintf('a2 = %.12g', C.a2)
            sprintf('b2 = %.12g', C.b2)
            sprintf('c2 = %.12g', C.c2)
            sprintf('delta = %.12g', delta_now)
            'delta = 5/6*b2^2 - a2*c2'
            'lattice minima:'
        }];

        for k = 1:numel(min_pos)
            lines{end+1} = sprintf('  min %d: psi2 = %.12g,  F = %.12g', ...
                k, min_pos(k), min_val(k));
        end

        if strcmpi(u_tag, 'mu') && isfinite(mu_folder_value)
            lines{end+1} = sprintf('mu(folder map) = %.12g', mu_folder_value);
        end

        coeff_box.Value = lines;
    end

    function render_empty_(msg)
        cla(axPsi);
        cla(axX);

        text(axPsi, 0.5, 0.5, msg, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 15);

        text(axX, 0.5, 0.5, msg, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 15);

        coeff_box.Value = {msg};
    end

    function apply_iqjq_()
        par_base.iq_pick = round(iq_input.Value);
        par_base.jq_pick = round(jq_input.Value);

        G2 = load_chi_grid_folderUI_(root, par_base.iq_pick, par_base.jq_pick);

        T_list  = G2.T_list;
        U_list  = G2.U_list;
        chi_map = G2.chi_map;
        dop_map = G2.doping_map;
        muf_map = G2.mu_folder_map;
        u_tag   = G2.u_tag;

        if isempty(T_list) || isempty(U_list)
            render_empty_('No valid chi grid after applying new iq/jq.');
            return;
        end

        [~, it_new] = min(abs(T_list - T_input.Value));
        iu_new = nearest_valid_in_row_(chi_map, it_new, 1);
        if ~isfinite(chi_map(it_new, iu_new))
            [~, iu_new] = find_first_valid_(chi_map);
        end

        U_slider.Limits = [1, max(1, numel(U_list))];
        U_slider.Value  = iu_new;

        refresh_labels_();
        update_();
    end
end

% =====================================================================
% loader
% =====================================================================
function G = load_chi_grid_folderUI_(root_dir, iq_pick, jq_pick)
    root_dir = string(root_dir);
    L = dir(fullfile(root_dir, "**", "chi*.txt"));

    iq_pick = round(iq_pick);
    jq_pick = round(jq_pick);

    n = numel(L);
    Tvals = nan(n,1);
    Uvals = nan(n,1);
    chiV  = nan(n,1);
    dopH  = nan(n,1);
    muF   = nan(n,1);

    u_tag = "";

    for k = 1:n
        fpath = string(fullfile(L(k).folder, L(k).name));

        [~,~,ukind,uval] = parse_TU_from_path_(fpath);
        if ukind=="" || ~isfinite(uval)
            continue;
        end
        if u_tag=="" && ukind~=""
            u_tag = ukind;
        end

        H = parse_header_Tdop_(fpath);
        if ~H.ok
            continue;
        end

        M = read_numeric_skiphash_(fpath);
        if isempty(M) || size(M,2) < 5
            continue;
        end

        Cc = detect_cols_(M);
        id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
        if isempty(id)
            continue;
        end

        Tvals(k) = round(H.T_K, 12);
        Uvals(k) = quantize_U_(ukind, uval);
        chiV(k)  = M(id, Cc.Re);
        dopH(k)  = round(H.doping, 12);

        if ukind=="mu"
            muF(k) = uval;
        end
    end

    if u_tag==""
        u_tag="doping";
    end

    mask = isfinite(Tvals) & isfinite(Uvals) & isfinite(chiV) & isfinite(dopH);
    Tvals = Tvals(mask);
    Uvals = Uvals(mask);
    chiV  = chiV(mask);
    dopH  = dopH(mask);
    muF   = muF(mask);

    T_list = sort(unique(Tvals));
    U_list = sort(unique(Uvals));

    NT = numel(T_list);
    NU = numel(U_list);

    chi_sum = zeros(NT,NU); chi_cnt = zeros(NT,NU);
    dop_sum = zeros(NT,NU); dop_cnt = zeros(NT,NU);
    mu_sum  = zeros(NT,NU); mu_cnt  = zeros(NT,NU);

    for i = 1:numel(Tvals)
        iT = find(T_list==Tvals(i), 1);
        iU = find(U_list==Uvals(i), 1);
        if isempty(iT) || isempty(iU)
            continue;
        end

        chi_sum(iT,iU) = chi_sum(iT,iU) + chiV(i);
        chi_cnt(iT,iU) = chi_cnt(iT,iU) + 1;

        dop_sum(iT,iU) = dop_sum(iT,iU) + dopH(i);
        dop_cnt(iT,iU) = dop_cnt(iT,iU) + 1;

        if isfinite(muF(i))
            mu_sum(iT,iU) = mu_sum(iT,iU) + muF(i);
            mu_cnt(iT,iU) = mu_cnt(iT,iU) + 1;
        end
    end

    chi_map = nan(NT,NU);
    m = chi_cnt > 0;
    chi_map(m) = chi_sum(m) ./ chi_cnt(m);

    doping_map = nan(NT,NU);
    md = dop_cnt > 0;
    doping_map(md) = dop_sum(md) ./ dop_cnt(md);

    mu_folder_map = nan(NT,NU);
    mm = mu_cnt > 0;
    mu_folder_map(mm) = mu_sum(mm) ./ mu_cnt(mm);

    G = struct( ...
        'T_list', T_list, ...
        'U_list', U_list, ...
        'chi_map', chi_map, ...
        'doping_map', doping_map, ...
        'mu_folder_map', mu_folder_map, ...
        'u_tag', u_tag);
end

function [Tfolder, Ufolder, ukind, uval] = parse_TU_from_path_(fpath)
    Tfolder = ""; Ufolder = ""; ukind = ""; uval = NaN; %#ok<NASGU>
    p = replace(string(fpath), "\", "/");

    tokMu = regexp(p, '/(mu[^/]+)/', 'tokens', 'once');
    if ~isempty(tokMu)
        Ufolder = string(tokMu{1});
        ukind = "mu";
        uval = extract_num_(Ufolder);
        if isfinite(uval)
            return;
        end
    end

    tokDp = regexp(p, '/(doping[^/]+)/', 'tokens', 'once');
    if ~isempty(tokDp)
        Ufolder = string(tokDp{1});
        ukind = "doping";
        uval = extract_num_(Ufolder);
    end
end

function x = extract_num_(s)
    tok = regexp(char(s), '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)', 'match', 'once');
    x = str2double(tok);
end

function u = quantize_U_(ukind, uval)
    if ukind=="doping"
        u = round(uval, 4);
    else
        u = round(uval, 12);
    end
end

function H = parse_header_Tdop_(fpath)
    H = struct('ok', false, 'T_K', NaN, 'doping', NaN);

    fid = fopen(fpath, 'r');
    if fid < 0
        return;
    end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';

    for t = 1:3000
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
            tok = regexp(s, ['^\#\s*T\s*=\s*' num], 'tokens', 'once');
            if isempty(tok)
                tok = regexp(s, ['^\#\s*T_K\s*=\s*' num], 'tokens', 'once');
            end
            if ~isempty(tok)
                H.T_K = str2double(tok{1});
            end
        end

        if ~isfinite(H.doping)
            tok = regexp(s, ['^\#\s*doping\s*=\s*' num], 'tokens', 'once');
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

function M = read_numeric_skiphash_(fpath)
    fid = fopen(fpath, 'r');
    if fid < 0
        M = [];
        return;
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

    lastFinite = find(any(isfinite(M),1), 1, 'last');
    if ~isempty(lastFinite)
        M = M(:,1:lastFinite);
    end
end

function C = detect_cols_(M)
    ncol = size(M,2);
    C = struct('iq', 1, 'jq', 2, 'Re', 5);

    if ncol >= 8
        C.iq = 2;
        C.jq = 3;
        C.Re = 6;
    end
end

function [it, iu] = find_first_valid_(M)
    [it, iu] = find(isfinite(M), 1, 'first');
    if isempty(it)
        it = 1;
        iu = 1;
    end
end

function iu2 = nearest_valid_in_row_(M, it, iu0)
    row = M(it,:);
    good = find(isfinite(row));
    if isempty(good)
        iu2 = iu0;
        return;
    end
    [~, k] = min(abs(good - iu0));
    iu2 = good(k);
end

function i = clamp_(i, N)
    i = max(1, min(N, i));
end

function imin = find_local_minima_(y)
    y = y(:);
    n = numel(y);
    imin = [];

    if n < 3
        return;
    end

    for i = 2:n-1
        if y(i) <= y(i-1) && y(i) <= y(i+1) && ...
           (y(i) < y(i-1) || y(i) < y(i+1))
            imin(end+1) = i; %#ok<AGROW>
        end
    end

    if y(1) < y(2)
        imin = [1, imin];
    end
    if y(end) < y(end-1)
        imin = [imin, n];
    end

    imin = unique(imin);
end