function H = hysteresis_width_vs_T_doping()
clc; close all;

par  = make_realchi_params(true);
coef = make_realchi_coeff(par);

% =========================
% USER SETTINGS
% =========================
T_list = [];                  % [] means auto-detect all T folders
sep_threshold = 1e-2;         % separation threshold for |psi_f|-|psi_b|
plot_each_T_hyst = false;

default_root = "E:/rg_master/data/";
if ~isfolder(default_root)
    default_root = string(pwd);
end

root = uigetdir(default_root, 'Select root folder that CONTAINS T folders');
if isequal(root,0), error('User cancelled.'); end
root = string(root);

Tfolders = detect_T_folders_(root);
if isempty(Tfolders)
    error("No T*.*** folders found directly under root: %s", root);
end

if isempty(T_list)
    T_list = [Tfolders.T].';
else
    T_list = T_list(:);
end

iq_pick = round(par.iq_pick);
jq_pick = round(par.jq_pick);

pmev = double(par.artificial_polar);
dmu  = 0.5 * pmev * 1e-3;

fprintf("Root: %s\n", root);
fprintf("Pick (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);
fprintf("artificial_polar=%.6g meV => dmu=%.6g eV\n", pmev, dmu);
fprintf("sep_threshold=%.3g\n", sep_threshold);

nT = numel(T_list);

Width_doping = nan(nT,1);
Dop_left     = nan(nT,1);
Dop_right    = nan(nT,1);
Max_sep      = nan(nT,1);
Res = cell(nT,1);

for it = 1:nT
    T_target = T_list(it);
    T_folder = find_T_folder_(Tfolders, T_target, par.T_tol);

    if strlength(T_folder) == 0
        warning("Skip T=%.6g K: no matched T folder.", T_target);
        continue;
    end

    fprintf("\n================ T = %.6g K ================\n", T_target);
    fprintf("T folder: %s\n", T_folder);

    try
        R = run_one_T_mu_simple_(T_folder, T_target, ...
            iq_pick, jq_pick, dmu, pmev, par, coef);
    catch ME
        warning("Skip T=%.6g K: %s", T_target, ME.message);
        continue;
    end

    mu_fwd  = R.mu_fwd(:);
    mu_bwd  = R.mu_bwd(:);
    dop_fwd = R.dop_fwd(:);
    dop_bwd = R.dop_bwd(:);

    psi_f = abs(R.psi_fwd(:,1));
    psi_b = abs(R.psi_bwd(:,1));

    [mu_b_sort, ordB] = sort(mu_bwd);
    psi_b_sort = psi_b(ordB);
    dop_b_sort = dop_bwd(ordB);

    psi_b_on_fwd = interp1(mu_b_sort, psi_b_sort, mu_fwd, 'linear', 'extrap');
    dop_b_on_fwd = interp1(mu_b_sort, dop_b_sort, mu_fwd, 'linear', 'extrap');

    sep = abs(psi_f - psi_b_on_fwd);
    Max_sep(it) = max(sep, [], 'omitnan');

    idx = find(sep > sep_threshold & isfinite(dop_fwd) & isfinite(dop_b_on_fwd));

    if isempty(idx)
        Width_doping(it) = 0;
        fprintf("[width] T=%.6g | no separation | max_sep=%.6g\n", ...
            T_target, Max_sep(it));
    else
        i1 = idx(1);
        i2 = idx(end);

        dop_sep = dop_fwd(idx);
        Dop_left(it)  = dop_sep(1);
        Dop_right(it) = dop_sep(end);

        Width_doping(it) = abs(Dop_right(it) - Dop_left(it));

        fprintf("[width] T=%.6g | dop_left=%.8g | dop_right=%.8g | width=%.8g | max_sep=%.6g\n", ...
            T_target, Dop_left(it), Dop_right(it), Width_doping(it), Max_sep(it));
    end

    R.sep_threshold = sep_threshold;
    R.sep_mu_grid = mu_fwd;
    R.sep_dop_grid = dop_fwd;
    R.sep_psi = sep;
    R.Width_doping = Width_doping(it);
    R.Dop_left = Dop_left(it);
    R.Dop_right = Dop_right(it);
    R.Max_sep = Max_sep(it);

    Res{it} = R;

    if plot_each_T_hyst
        plot_one_T_hyst_with_sep_(R, par);
    end
end

fs = par.plot.fontSize;

figure('Color','w','Units','pixels','Position',[120 120 760 560]);
plot(T_list, Width_doping, 'o-', 'LineWidth',2, 'MarkerSize',7);
grid on; box on;
xlabel('$T$ (K)', 'Interpreter','latex');
ylabel('$\Delta n_{\rm hyst}$ $(10^{12}\,\mathrm{cm}^{-2})$', 'Interpreter','latex');
title(sprintf('Hysteresis width from $|\\Delta\\psi|>%.1e$', sep_threshold), ...
    'Interpreter','latex','FontWeight','normal');
set(gca,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1.2,'TickDir','out');

figure('Color','w','Units','pixels','Position',[160 160 760 560]);
plot(T_list, Max_sep, 'o-', 'LineWidth',2, 'MarkerSize',7);
grid on; box on;
xlabel('$T$ (K)', 'Interpreter','latex');
ylabel('$\max |\Delta\psi|$', 'Interpreter','latex');
title('Maximum forward-backward separation', ...
    'Interpreter','latex','FontWeight','normal');
set(gca,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1.2,'TickDir','out');

H = struct();
H.root = root;
H.T_list = T_list;
H.sep_threshold = sep_threshold;
H.Width_doping = Width_doping;
H.Dop_left = Dop_left;
H.Dop_right = Dop_right;
H.Max_sep = Max_sep;
H.results = Res;
H.iq_pick = iq_pick;
H.jq_pick = jq_pick;
H.artificial_polar = pmev;
H.dmu_eV = dmu;
H.par = par;

end

function Tfolders = detect_T_folders_(root)

D = dir(root);
D = D([D.isdir]);

Ts = [];
Paths = strings(0,1);

for i = 1:numel(D)
    name = string(D(i).name);
    if name == "." || name == ".."
        continue;
    end

    tok = regexp(name, "^T([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)$", "tokens", "once");
    if isempty(tok)
        continue;
    end

    Tval = str2double(tok{1});
    if ~isfinite(Tval)
        continue;
    end

    Ts(end+1,1) = Tval; %#ok<AGROW>
    Paths(end+1,1) = string(fullfile(root, name)); %#ok<AGROW>
end

[Ts, ord] = sort(Ts);
Paths = Paths(ord);

Tfolders = struct('T', num2cell(Ts), 'path', cellstr(Paths));

end

function T_folder = find_T_folder_(Tfolders, T_target, T_tol)

Tall = [Tfolders.T].';
[diffMin, id] = min(abs(Tall - T_target));

if isempty(id) || diffMin > max(T_tol, 1e-8)
    T_folder = "";
else
    T_folder = string(Tfolders(id).path);
end

end

function R = run_one_T_mu_simple_(T_folder, T_target, iq_pick, jq_pick, dmu, pmev, par, coef)

B = collect_mu_dop_chi_curve_(T_folder, iq_pick, jq_pick);

mu0  = B.mu(:);
dop0 = B.dop(:);
chi0 = B.chi(:);

fprintf("[baseline] points=%d | mu=[%.6g, %.6g] | dop=[%.6g, %.6g]\n", ...
    numel(mu0), min(mu0), max(mu0), min(dop0), max(dop0));

if numel(mu0) < 5
    error("Too few points in this T folder.");
end

[mu0, chi0, dop0] = collapse_same_mu_(mu0, chi0, dop0);

Fchi0 = griddedInterpolant(mu0, chi0, "linear", "none");
Fdop0 = griddedInterpolant(mu0, dop0, "linear", "none");

if abs(dmu) < 1e-18
    chi_art_mu = chi0;
    ok = true(size(mu0));
else
    chi_art_mu = 0.5 * (Fchi0(mu0 + dmu) + Fchi0(mu0 - dmu));
    ok = isfinite(chi_art_mu);
end

mu_use  = mu0(ok);
chi_use = chi_art_mu(ok);
dop_use = Fdop0(mu_use);

ok2 = isfinite(mu_use) & isfinite(chi_use) & isfinite(dop_use);
mu_use  = mu_use(ok2);
chi_use = chi_use(ok2);
dop_use = dop_use(ok2);

if numel(mu_use) < 5
    error("After shift-average, too few finite points.");
end

[mu_use, ord] = sort(mu_use);
chi_use = chi_use(ord);
dop_use = dop_use(ord);

[mu_use, chi_use, dop_use] = collapse_same_mu_(mu_use, chi_use, dop_use);

Fchi_art = griddedInterpolant(mu_use, chi_use, "linear", "none");
Fdop_art = griddedInterpolant(mu_use, dop_use, "linear", "none");

mu_fwd = mu_use(:);
mu_bwd = flipud(mu_use(:));
N = numel(mu_fwd);

fprintf("[mu scan] forward: [%.6g, %.6g], backward: [%.6g, %.6g], N=%d\n", ...
    mu_fwd(1), mu_fwd(end), mu_bwd(1), mu_bwd(end), N);

psi_f = nan(N,2);
psi_b = nan(N,2);
a1_f  = nan(N,1);
a1_b  = nan(N,1);
dop_f = nan(N,1);
dop_b = nan(N,1);
chi_f = nan(N,1);
chi_b = nan(N,1);

prev = [par.hyst.psi1_0, par.hyst.psi2_0];

for i = 1:N
    mu = mu_fwd(i);
    chi_used = Fchi_art(mu);
    dop = Fdop_art(mu);

    C = coef.eval(T_target, dop, chi_used);

    [p1,p2,~] = minimize_free_energy( ...
        C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda, ...
        prev(1), prev(2), false);

    pick = [p1,p2];
    psi_f(i,:) = pick;
    prev = pick;

    a1_f(i)  = C.a1;
    dop_f(i) = dop;
    chi_f(i) = chi_used;
end

prev = psi_f(end,:);

for i = 1:N
    mu = mu_bwd(i);
    chi_used = Fchi_art(mu);
    dop = Fdop_art(mu);

    C = coef.eval(T_target, dop, chi_used);

    [p1,p2,~] = minimize_free_energy( ...
        C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda, ...
        prev(1), prev(2), false);

    pick = [p1,p2];
    psi_b(i,:) = pick;
    prev = pick;

    a1_b(i)  = C.a1;
    dop_b(i) = dop;
    chi_b(i) = chi_used;
end

R = struct();
R.T = T_target;
R.T_folder = T_folder;
R.iq_pick = iq_pick;
R.jq_pick = jq_pick;
R.artificial_polar = pmev;
R.dmu_eV = dmu;

R.mu_raw  = mu0;
R.dop_raw = dop0;
R.chi_raw = chi0;

R.mu_used  = mu_use;
R.dop_used = dop_use;
R.chi_used = chi_use;

R.mu_fwd = mu_fwd;
R.mu_bwd = mu_bwd;

R.dop_fwd = dop_f;
R.dop_bwd = dop_b;

R.psi_fwd = psi_f;
R.psi_bwd = psi_b;

R.a1_fwd = a1_f;
R.a1_bwd = a1_b;

R.chi_fwd = chi_f;
R.chi_bwd = chi_b;

end

function B = collect_mu_dop_chi_curve_(T_folder, iq_pick, jq_pick)

L = dir(fullfile(T_folder, "**", "chi*.txt"));
L = L(~[L.isdir]);

if isempty(L)
    error("No chi*.txt under T folder: %s", T_folder);
end

mu  = nan(numel(L),1);
dop = nan(numel(L),1);
chi = nan(numel(L),1);

n = 0;
dropH = 0;
dropTbl = 0;
dropQ = 0;

for k = 1:numel(L)
    fpath = string(fullfile(L(k).folder, L(k).name));
    folder = string(L(k).folder);

    H = parse_header_mu_dop_(fpath);

    mu_path = parse_mu_from_path_(folder);
    mu_use = mu_path;

    if ~isfinite(mu_use)
        mu_use = H.mu;
    end

    if ~isfinite(mu_use) || ~isfinite(H.dop)
        dropH = dropH + 1;
        continue;
    end

    try
        M = readmatrix(fpath, "FileType","text", "CommentStyle","#");
    catch
        dropTbl = dropTbl + 1;
        continue;
    end

    if isempty(M) || size(M,2) < 6
        dropTbl = dropTbl + 1;
        continue;
    end

    Cc = detect_cols_(M);
    id = find(M(:,Cc.iq) == iq_pick & M(:,Cc.jq) == jq_pick, 1);

    if isempty(id)
        dropQ = dropQ + 1;
        continue;
    end

    n = n + 1;
    mu(n)  = mu_use;
    dop(n) = H.dop;
    chi(n) = M(id, Cc.Re);
end

mu  = mu(1:n);
dop = dop(1:n);
chi = chi(1:n);

ok = isfinite(mu) & isfinite(dop) & isfinite(chi);
mu  = mu(ok);
dop = dop(ok);
chi = chi(ok);

fprintf("[collect] files=%d kept=%d | dropH=%d dropTbl=%d dropQ=%d\n", ...
    numel(L), numel(mu), dropH, dropTbl, dropQ);

B = struct('mu',mu,'dop',dop,'chi',chi);

end

function H = parse_header_mu_dop_(fpath)

H = struct('mu',NaN,'dop',NaN);

fid = fopen(fpath,'r');
if fid < 0
    return;
end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

for t = 1:3000
    ln = fgetl(fid);
    if ~ischar(ln), break; end

    s0 = strtrim(string(ln));
    if ~startsWith(s0,"#"), break; end

    s = strtrim(erase(s0,"#"));

    if ~isfinite(H.dop)
        tok = regexp(s, "(?:^|\s)(?:doping|dop)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.dop = str2double(tok{1}); end
    end

    if ~isfinite(H.mu)
        tok = regexp(s, "(?:^|\s)(?:mu|EF|E_F)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.mu = str2double(tok{1}); end
    end
end

end

function mu = parse_mu_from_path_(folder)

folder = replace(string(folder),"\","/");
num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

tok = regexp(folder, "/mu\s*=?\s*" + num + "(?:/|$)", "tokens","once","ignorecase");

if ~isempty(tok)
    mu = str2double(tok{1});
else
    mu = NaN;
end

end

function C = detect_cols_(M)

if size(M,2) >= 8
    C = struct('iq',2,'jq',3,'Re',6);
else
    C = struct('iq',1,'jq',2,'Re',5);
end

end

function [mu2, chi2, dop2] = collapse_same_mu_(mu, chi, dop)

mu  = double(mu(:));
chi = double(chi(:));
dop = double(dop(:));

ok = isfinite(mu) & isfinite(chi) & isfinite(dop);
mu  = mu(ok);
chi = chi(ok);
dop = dop(ok);

[mu_u,~,ic] = unique(mu);
chi_m = accumarray(ic, chi, [], @mean);
dop_m = accumarray(ic, dop, [], @mean);

[mu2,ord] = sort(mu_u);
chi2 = chi_m(ord);
dop2 = dop_m(ord);

end

function plot_one_T_hyst_with_sep_(R, par)

fs = par.plot.fontSize;

figure('Color','w','Units','pixels','Position',[80 80 720 560]);
plot(R.dop_fwd, abs(R.psi_fwd(:,1)), '-',  'LineWidth',2); hold on;
plot(R.dop_bwd, abs(R.psi_bwd(:,1)), '--', 'LineWidth',2);

if isfinite(R.Dop_left) && isfinite(R.Dop_right)
    xline(R.Dop_left,  ':', 'LineWidth',1.5);
    xline(R.Dop_right, ':', 'LineWidth',1.5);
end

grid on; box on;
xlabel('doping $(10^{12}\,\mathrm{cm}^{-2})$','Interpreter','latex');
ylabel('$|\psi|$','Interpreter','latex');
title(sprintf('$|\\psi|$ hysteresis, $T=%.6g$ K, width $=%.4g$', ...
    R.T, R.Width_doping), 'Interpreter','latex','FontWeight','normal');
legend({'forward','backward'}, 'Interpreter','latex', 'Location','best');
set(gca,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1.2,'TickDir','out');

end