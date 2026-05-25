%% plot_free_energy_manual_mu.m
% 直接运行脚本：
%   - 手动在脚本里指定 root / T_target / mu_list
%   - 自动读取 chi 数据
%   - 自动找最接近的 (T, mu)
%   - 画两个 figure:
%       1) 所有 mu 的 electronic free energy F_psi(psi1)
%       2) 所有 mu 的 lattice free energy    F_X(psi2)

clear; clc;

%% =========================================================
% 1. parameters
% ==========================================================
par = make_realchi_params(true);

% ---- default root ----
default_root = "E:/rg_master/data";
if isfield(par,'io') && isfield(par.io,'default_root') ...
        && isfolder(par.io.default_root)
    default_root = string(par.io.default_root);
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

T_target = 6.5;
mu_list  = [0.58, 0.7, 0.7885];   % <<<<<< 在这里设置你想画的 mu 列表

% 如果你想覆盖默认 q 点，就改这里；否则就用 make_realchi_params 里的默认值
% par.iq_pick = -267;
% par.jq_pick = -267;

% 如果你想手动改 alpha，也可以直接改
% par.Lat_alpha = 0.2;

fprintf('Root       = %s\n', root);
fprintf('T_target   = %.12g\n', T_target);
fprintf('mu_list    = [ ');
fprintf('%.12g ', mu_list);
fprintf(']\n');
fprintf('iq_pick    = %d\n', par.iq_pick);
fprintf('jq_pick    = %d\n', par.jq_pick);
fprintf('Lat_alpha  = %.12g\n', par.Lat_alpha);

%% =========================================================
% 3. load chi grid
% ==========================================================
G = load_chi_grid_folderUI_(root, par.iq_pick, par.jq_pick);

T_list  = G.T_list;
U_list  = G.U_list;
chi_map = G.chi_map;
dop_map = G.doping_map;
muf_map = G.mu_folder_map;
u_tag   = G.u_tag;

if isempty(T_list) || isempty(U_list)
    error('No valid chi grid found under selected folder.');
end

if ~strcmpi(u_tag, 'mu')
    error('This script requires mu-folder dataset, but detected type is "%s".', u_tag);
end

%% =========================================================
% 4. fixed x-ranges
% ==========================================================
if isfield(par, 'psi1_lim') && numel(par.psi1_lim)==2
    x1 = linspace(par.psi1_lim(1), par.psi1_lim(2), 2001);
    xlim1 = par.psi1_lim;
else
    x1 = linspace(-10, 10, 2001);
    xlim1 = [x1(1), x1(end)];
end

if isfield(par, 'psi2_lim') && numel(par.psi2_lim)==2
    x2 = linspace(par.psi2_lim(1), par.psi2_lim(2), 4001);
    xlim2 = par.psi2_lim;
else
    x2 = linspace(-20, 20, 4001);
    xlim2 = [x2(1), x2(end)];
end

%% =========================================================
% 5. find nearest T
% ==========================================================
[~, it] = min(abs(T_list - T_target));
T_now = T_list(it);

coef = make_realchi_coeff(par);

% 存结果，方便后面打印
res = struct('mu_req', {}, 'mu_used', {}, 'mu_folder_map', {}, ...
             'doping', {}, 'chi', {}, 'C', {}, 'delta', {}, ...
             'y1', {}, 'y2', {}, 'idx1', {}, ...
             'min_pos2', {}, 'min_val2', {});

%% =========================================================
% 6. loop over mu_list
% ==========================================================
for im = 1:numel(mu_list)
    mu_target = mu_list(im);

    [~, iu] = min(abs(U_list - mu_target));

    if ~isfinite(chi_map(it, iu))
        iu = nearest_valid_in_row_(chi_map, it, iu);
    end

    if ~isfinite(chi_map(it, iu))
        warning('No valid point near T = %.12g and mu = %.12g. Skip.', T_target, mu_target);
        continue;
    end

    mu_now        = U_list(iu);
    chi_used      = chi_map(it, iu);
    dop_eval      = dop_map(it, iu);
    mu_folder_val = muf_map(it, iu);

    C = coef.eval(T_now, dop_eval, chi_used);
    delta_now = (5/6) * C.b2^2 - C.a2 * C.c2;

    % free energies
    Fpsi = @(psi) 0.5 * C.a1 * psi.^2 + (1/factorial(4)) * C.b1 * psi.^4;
    FX   = @(X)   0.5 * C.a2 * X.^2 ...
                + (1/factorial(4)) * C.b2 * X.^4 ...
                + (1/factorial(6)) * C.c2 * X.^6;

    y1 = Fpsi(x1);
    y2 = FX(x2);

    [~, idx1] = min(y1);

    imin2 = find_local_minima_(y2);
    if isempty(imin2)
        [~, idx2] = min(y2);
        imin2 = idx2;
    end

    min_pos2 = x2(imin2);
    min_val2 = y2(imin2);
    [min_pos2, order2] = sort(min_pos2);
    min_val2 = min_val2(order2);

    res(end+1).mu_req = mu_target; %#ok<SAGROW>
    res(end).mu_used = mu_now;
    res(end).mu_folder_map = mu_folder_val;
    res(end).doping = dop_eval;
    res(end).chi = chi_used;
    res(end).C = C;
    res(end).delta = delta_now;
    res(end).y1 = y1;
    res(end).y2 = y2;
    res(end).idx1 = idx1;
    res(end).min_pos2 = min_pos2;
    res(end).min_val2 = min_val2;
end

if isempty(res)
    error('No valid mu points found for the requested mu_list at T = %.12g.', T_target);
end

%% =========================================================
% 7. figure 1: all electronic curves
% ==========================================================
figure('Name', 'Electronic free energy (all mu)', 'Color', 'w');
hold on;

leg1 = cell(1, numel(res));
for k = 1:numel(res)
    plot(x1, res(k).y1, 'LineWidth', 2);
    plot(x1(res(k).idx1), res(k).y1(res(k).idx1), 'o', ...
        'MarkerSize', 6, 'LineWidth', 1.2, 'MarkerFaceColor', 'none');

    leg1{k} = sprintf('\\mu_{req}=%.4g, \\mu_{used}=%.4g', ...
        res(k).mu_req, res(k).mu_used);
end

hold off;
grid on;
box on;
xlim(xlim1);
xlabel('$\psi_1$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$F_{\psi}(\psi_1)$', 'Interpreter', 'latex', 'FontSize', 16);
title(sprintf('$F_{\\psi}$ at $T=%.6g$', T_now), ...
    'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
legend(leg1, 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', 14);

%% =========================================================
% 8. figure 2: all lattice curves
% ==========================================================
figure('Name', 'Lattice free energy (all mu)', 'Color', 'w');
hold on;

leg2 = cell(1, numel(res));
for k = 1:numel(res)
    plot(x2, res(k).y2, 'LineWidth', 2);
    plot(res(k).min_pos2, res(k).min_val2, 'o', ...
        'MarkerSize', 6, 'LineWidth', 1.2, 'MarkerFaceColor', 'none');

    leg2{k} = sprintf('\\mu_{req}=%.4g, \\mu_{used}=%.4g', ...
        res(k).mu_req, res(k).mu_used);
end

hold off;
grid on;
box on;
xlim(xlim2);
xlabel('$\psi_2$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$F_X(\psi_2)$', 'Interpreter', 'latex', 'FontSize', 16);
title(sprintf('$F_X$ at $T=%.6g$', T_now), ...
    'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
legend(leg2, 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', 14);

%% =========================================================
% 9. summary print
% ==========================================================
fprintf('\n================ free energy summary ================\n');
fprintf('iq = %d, jq = %d\n', par.iq_pick, par.jq_pick);
fprintf('T requested = %.12g\n', T_target);
fprintf('T used      = %.12g\n', T_now);
fprintf('alpha       = %.12g\n', par.Lat_alpha);

for k = 1:numel(res)
    fprintf('\n---- mu case %d ----\n', k);
    fprintf('mu requested     = %.12g\n', res(k).mu_req);
    fprintf('mu used          = %.12g\n', res(k).mu_used);
    fprintf('mu(folder map)   = %.12g\n', res(k).mu_folder_map);
    fprintf('doping(header)   = %.12g\n', res(k).doping);
    fprintf('chi              = %.12g\n', res(k).chi);
    fprintf('a1               = %.12g\n', res(k).C.a1);
    fprintf('b1               = %.12g\n', res(k).C.b1);
    fprintf('a2               = %.12g\n', res(k).C.a2);
    fprintf('b2               = %.12g\n', res(k).C.b2);
    fprintf('c2               = %.12g\n', res(k).C.c2);
    fprintf('delta            = %.12g\n', res(k).delta);
    for j = 1:numel(res(k).min_pos2)
        fprintf('lattice min %d: psi2 = %.12g, F = %.12g\n', ...
            j, res(k).min_pos2(j), res(k).min_val2(j));
    end
end
fprintf('=====================================================\n');

%% =====================================================================
% local functions
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