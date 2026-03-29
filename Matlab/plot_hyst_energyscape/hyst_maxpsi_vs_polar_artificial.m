function OUT = hyst_maxpsi_vs_polar_cached()
% hyst_maxpsi_vs_polar_cached
% ------------------------------------------------------------
% Read chi data ONLY ONCE, then scan custom artificial polar list.
% For each polar:
%   1) build chi_art(mu) = 0.5*(chi(mu+dmu)+chi(mu-dmu))
%   2) map to doping
%   3) run forward/backward hysteresis continuation
%   4) extract max |psi1|
%
% Depends on existing core functions:
%   - make_realchi_params
%   - make_realchi_coeff
%   - free_energy
%   - minimize_free_energy
%
% No modification to your original hysteresis script is needed.
% ------------------------------------------------------------

% ---------------- params / coeff ----------------
par  = make_realchi_params(true);
coef = make_realchi_coeff(par);

% ---------------- choose folder -----------------
default_root = "/Users/haoranyan/rg_master/data/";
if ~isfolder(default_root)
    default_root = string(pwd);
end
root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
if isequal(root,0), error('User cancelled.'); end
root = string(root);

fprintf("Root: %s\n", root);

% ---------------- custom polar list -------------
polarList_meV = 0:0.5:8;   % <<< 改这里

% ---------------- q point / T -------------------
iq_pick  = round(par.iq_pick);
jq_pick  = round(par.jq_pick);
T_target = par.T_target;
T_tol    = par.T_tol;

fprintf("Pick (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);
fprintf("T_target=%.6g K, T_tol=%.3g K\n", T_target, T_tol);

% =========================================================
% 1) collect baseline curve once
% =========================================================
B = collect_mu_dop_chi_curve_(root, T_target, T_tol, iq_pick, jq_pick);

mu0  = B.mu(:);
dop0 = B.dop(:);
chi0 = B.chi(:);

fprintf("[baseline] points=%d | mu=[%.6g, %.6g] | dop=[%.6g, %.6g]\n", ...
    numel(mu0), min(mu0), max(mu0), min(dop0), max(dop0));

if numel(mu0) < 5
    error("Too few points at this T. Check T_tol / header parse / (iq,jq).");
end

[mu0, chi0, dop0] = collapse_same_mu_(mu0, chi0, dop0);

Fchi0 = griddedInterpolant(mu0, chi0, "linear", "none");
Fdop0 = griddedInterpolant(mu0, dop0, "linear", "none");

% ---------------- build common doping path once ----------------
[dop_fwd, dop_bwd] = build_hyst_doping_grid_(dop0, par.hyst);
Nscan = numel(dop_fwd);

fprintf("[scan window] forward: [%.6g, %.6g], backward: [%.6g, %.6g], N=%d\n", ...
    dop_fwd(1), dop_fwd(end), dop_bwd(1), dop_bwd(end), Nscan);

% =========================================================
% 2) loop over polar
% =========================================================
np = numel(polarList_meV);

maxPsi_fwd = nan(np,1);
maxPsi_bwd = nan(np,1);
maxPsi_all = nan(np,1);
Hall = cell(np,1);

for ip = 1:np
    pmev = double(polarList_meV(ip));
    dmu  = 0.5 * pmev * 1e-3;   % meV -> eV, then /2

    fprintf("\n==================================================\n");
    fprintf("[polar %d/%d] artificial_polar = %.6g meV  => dmu = %.6g eV\n", ...
        ip, np, pmev, dmu);

    % -----------------------------------------------------
    % build chi_art(mu), then map to doping
    % -----------------------------------------------------
    if abs(dmu) < 1e-18
        mu_use  = mu0;
        chi_use = chi0;
        dop_use = dop0;
    else
        chi_art_mu = 0.5 * (Fchi0(mu0 + dmu) + Fchi0(mu0 - dmu));
        dop_art_mu = Fdop0(mu0);

        ok = isfinite(mu0) & isfinite(chi_art_mu) & isfinite(dop_art_mu);
        mu_use  = mu0(ok);
        chi_use = chi_art_mu(ok);
        dop_use = dop_art_mu(ok);
    end

    if numel(mu_use) < 5
        warning("polar=%.6g meV skipped: too few finite shifted points.", pmev);
        continue;
    end

    [dop_use, ord] = sort(dop_use);
    chi_use = chi_use(ord);
    mu_use  = mu_use(ord);

    [dop_use, chi_use, mu_use] = collapse_same_dop_(dop_use, chi_use, mu_use);

    if numel(dop_use) < 5
        warning("polar=%.6g meV skipped: too few unique doping points.", pmev);
        continue;
    end

    chi_of_dop = @(d) interp1(dop_use, chi_use, d, 'linear', 'extrap');

    % -----------------------------------------------------
    % forward continuation
    % -----------------------------------------------------
    psi_fwd = nan(Nscan,2);
    a1_fwd  = nan(Nscan,1);
    chi_fwd = nan(Nscan,1);

    prev = [par.hyst.psi1_0, par.hyst.psi2_0];

    for i = 1:Nscan
        dop = dop_fwd(i);
        chi_now = chi_of_dop(dop);
        C = coef.eval(T_target, dop, chi_now);

        [p1,p2,~] = minimize_free_energy( ...
            C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda, ...
            prev(1), prev(2), false);

        psi_fwd(i,:) = [p1,p2];
        prev = [p1,p2];

        a1_fwd(i)  = C.a1;
        chi_fwd(i) = chi_now;
    end

    % -----------------------------------------------------
    % backward continuation
    % -----------------------------------------------------
    psi_bwd = nan(Nscan,2);
    a1_bwd  = nan(Nscan,1);
    chi_bwd = nan(Nscan,1);

    prev = psi_fwd(end,:);

    for i = 1:Nscan
        dop = dop_bwd(i);
        chi_now = chi_of_dop(dop);
        C = coef.eval(T_target, dop, chi_now);

        [p1,p2,~] = minimize_free_energy( ...
            C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda, ...
            prev(1), prev(2), false);

        psi_bwd(i,:) = [p1,p2];
        prev = [p1,p2];

        a1_bwd(i)  = C.a1;
        chi_bwd(i) = chi_now;
    end

    % -----------------------------------------------------
    % optional plotting-only Lorentz smoothing
    % -----------------------------------------------------
    if isfield(par, 'smooth')
        [psi1_f_plot, psi1_b_plot, gamma_used] = apply_lorentz_on_doping_( ...
            dop_fwd, psi_fwd(:,1), dop_bwd, psi_bwd(:,1), par.smooth);
    else
        psi1_f_plot = psi_fwd(:,1);
        psi1_b_plot = psi_bwd(:,1);
        gamma_used = NaN;
    end

    % -----------------------------------------------------
    % max electronic order
    % -----------------------------------------------------
    maxPsi_fwd(ip) = max(abs(psi_fwd(:,1)), [], 'omitnan');
    maxPsi_bwd(ip) = max(abs(psi_bwd(:,1)), [], 'omitnan');
    maxPsi_all(ip) = max([maxPsi_fwd(ip), maxPsi_bwd(ip)], [], 'omitnan');

    fprintf("max |psi1| forward  = %.6g\n", maxPsi_fwd(ip));
    fprintf("max |psi1| backward = %.6g\n", maxPsi_bwd(ip));
    fprintf("max |psi1| overall  = %.6g\n", maxPsi_all(ip));

    H = struct();
    H.polar_meV  = pmev;
    H.dmu_eV     = dmu;
    H.mu_used    = mu_use;
    H.dop_used   = dop_use;
    H.chi_used   = chi_use;
    H.dop_fwd    = dop_fwd;
    H.dop_bwd    = dop_bwd;
    H.psi_fwd    = psi_fwd;
    H.psi_bwd    = psi_bwd;
    H.psi1_f_plot = psi1_f_plot;
    H.psi1_b_plot = psi1_b_plot;
    H.gamma_used = gamma_used;
    H.a1_fwd     = a1_fwd;
    H.a1_bwd     = a1_bwd;
    H.chi_fwd    = chi_fwd;
    H.chi_bwd    = chi_bwd;
    H.maxPsi_fwd = maxPsi_fwd(ip);
    H.maxPsi_bwd = maxPsi_bwd(ip);
    H.maxPsi_all = maxPsi_all(ip);

    Hall{ip} = H;
end

% =========================================================
% 3) plot result
% =========================================================
fig = figure('Color','w','Position',[100 100 850 620]);
plot(polarList_meV, maxPsi_fwd, 'o-','LineWidth',2,'MarkerSize',7); hold on;
plot(polarList_meV, maxPsi_bwd, 's--','LineWidth',2,'MarkerSize',7);
plot(polarList_meV, maxPsi_all, 'd-','LineWidth',2.2,'MarkerSize',7);
hold off;
grid on; box on;

xlabel('artificial polar (meV)');
ylabel('max |\psi_1|');
title(sprintf('Max electronic order vs polar, T = %.6g K, q = (%d,%d)', ...
    T_target, iq_pick, jq_pick));
legend({'forward','backward','overall'}, 'Location','best');
set(gca,'FontSize',16,'LineWidth',1);

% =========================================================
% 4) save
% =========================================================
out_dir = fullfile(root, "..", "plot");
if ~exist(out_dir, "dir")
    mkdir(out_dir);
end

out_png = fullfile(out_dir, sprintf( ...
    "maxPsi_vs_polar_T%.4f_iq%d_jq%d.png", ...
    T_target, iq_pick, jq_pick));

out_mat = fullfile(out_dir, sprintf( ...
    "maxPsi_vs_polar_T%.4f_iq%d_jq%d.mat", ...
    T_target, iq_pick, jq_pick));

exportgraphics(fig, out_png, 'Resolution', 300);
save(out_mat, 'polarList_meV', 'maxPsi_fwd', 'maxPsi_bwd', 'maxPsi_all', ...
    'Hall', 'root', 'T_target', 'iq_pick', 'jq_pick', 'par');

fprintf("\nSaved figure: %s\n", out_png);
fprintf("Saved data  : %s\n", out_mat);

OUT = struct();
OUT.polarList_meV = polarList_meV(:);
OUT.maxPsi_fwd = maxPsi_fwd;
OUT.maxPsi_bwd = maxPsi_bwd;
OUT.maxPsi_all = maxPsi_all;
OUT.Hall = Hall;
OUT.png = out_png;
OUT.mat = out_mat;

end

% =====================================================================
% collect chi data once
% =====================================================================
function B = collect_mu_dop_chi_curve_(root, T_target, T_tol, iq_pick, jq_pick)

L = dir(fullfile(root, "**", "chi*.txt"));
L = L(~[L.isdir]);

if isempty(L)
    error("No chi*.txt found under: %s", root);
end

Tlo = T_target - T_tol;
Thi = T_target + T_tol;

mu  = nan(numel(L),1);
dop = nan(numel(L),1);
chi = nan(numel(L),1);

n = 0;
dropT = 0; dropH = 0; dropTbl = 0; dropQ = 0;

for k = 1:numel(L)
    fpath = string(fullfile(L(k).folder, L(k).name));
    folder = string(L(k).folder);

    T_path = parse_T_from_path_(folder);
    H = parse_header_mu_dop_T_(fpath);

    if isfinite(T_path)
        T_use = T_path;
    else
        T_use = H.T;
    end

    if ~isfinite(T_use)
        dropH = dropH + 1;
        continue;
    end

    if ~(T_use >= Tlo && T_use <= Thi)
        dropT = dropT + 1;
        continue;
    end

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

    C = detect_cols_(M);
    id = find(M(:,C.iq)==iq_pick & M(:,C.jq)==jq_pick, 1);

    if isempty(id)
        dropQ = dropQ + 1;
        continue;
    end

    n = n + 1;
    mu(n)  = mu_use;
    dop(n) = H.dop;
    chi(n) = M(id, C.Re);
end

mu  = mu(1:n);
dop = dop(1:n);
chi = chi(1:n);

ok = isfinite(mu) & isfinite(dop) & isfinite(chi);
mu = mu(ok);
dop = dop(ok);
chi = chi(ok);

fprintf("[collect] files=%d kept=%d | dropT=%d dropH=%d dropTbl=%d dropQ=%d\n", ...
    numel(L), numel(mu), dropT, dropH, dropTbl, dropQ);

B = struct('mu',mu,'dop',dop,'chi',chi);

end

function H = parse_header_mu_dop_T_(fpath)

H = struct('T',NaN,'mu',NaN,'dop',NaN);

fid = fopen(fpath,'r');
if fid < 0
    return;
end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';

for t = 1:3000
    ln = fgetl(fid);
    if ~ischar(ln), break; end

    s0 = strtrim(string(ln));
    if ~startsWith(s0, "#"), break; end
    s = strtrim(erase(s0, "#"));

    if ~isfinite(H.T)
        tok = regexp(s, "(?:^|\s)T(?:_K)?\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.T = str2double(tok{1}); end
    end

    if ~isfinite(H.dop)
        tok = regexp(s, "(?:^|\s)(?:doping|dop)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.dop = str2double(tok{1}); end
    end

    if ~isfinite(H.mu)
        tok = regexp(s, "(?:^|\s)(?:mu|mu_eV|EF|EF_eV|E_F)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.mu = str2double(tok{1}); end
    end
end

end

function T = parse_T_from_path_(folder)
folder = replace(string(folder), "\", "/");
num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';

tok = regexp(folder, "/T\s*=?\s*" + num + "(?:/|$)", "tokens","once","ignorecase");
if ~isempty(tok)
    T = str2double(tok{1});
    return;
end

tok = regexp(folder, "/temp\s*=?\s*" + num + "(?:/|$)", "tokens","once","ignorecase");
if ~isempty(tok)
    T = str2double(tok{1});
else
    T = NaN;
end
end

function mu = parse_mu_from_path_(folder)
folder = replace(string(folder), "\", "/");
num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';

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
mu = mu(ok);
chi = chi(ok);
dop = dop(ok);

[mu_u,~,ic] = unique(mu);
chi_m = accumarray(ic, chi, [], @mean);
dop_m = accumarray(ic, dop, [], @mean);

[mu2,ord] = sort(mu_u);
chi2 = chi_m(ord);
dop2 = dop_m(ord);
end

function [dop2, chi2, mu2] = collapse_same_dop_(dop, chi, mu)
dop = double(dop(:));
chi = double(chi(:));
mu  = double(mu(:));

ok = isfinite(dop) & isfinite(chi) & isfinite(mu);
dop = dop(ok);
chi = chi(ok);
mu  = mu(ok);

[dop_u,~,ic] = unique(dop);
chi_m = accumarray(ic, chi, [], @mean);
mu_m  = accumarray(ic, mu,  [], @mean);

[dop2,ord] = sort(dop_u);
chi2 = chi_m(ord);
mu2  = mu_m(ord);
end

function [dop_fwd, dop_bwd] = build_hyst_doping_grid_(dop0, hyst)
dop0 = double(dop0(:));
dop0 = dop0(isfinite(dop0));

if isempty(dop0)
    error("build_hyst_doping_grid_: dop0 is empty.");
end

dmin_all = min(dop0);
dmax_all = max(dop0);

mode = "N";
if isfield(hyst,"grid_mode") && ~isempty(hyst.grid_mode)
    mode = string(hyst.grid_mode);
end

direction = "ascend";
if isfield(hyst,"forward_direction") && ~isempty(hyst.forward_direction)
    direction = lower(string(hyst.forward_direction));
end

scan_start = NaN;
scan_end   = NaN;

if isfield(hyst,"scan_start") && ~isempty(hyst.scan_start)
    scan_start = double(hyst.scan_start);
end
if isfield(hyst,"scan_end") && ~isempty(hyst.scan_end)
    scan_end = double(hyst.scan_end);
end

if ~isfinite(scan_start), scan_start = dmin_all; end
if ~isfinite(scan_end),   scan_end   = dmax_all; end

dlo_req = min(scan_start, scan_end);
dhi_req = max(scan_start, scan_end);

dlo = max(dmin_all, dlo_req);
dhi = min(dmax_all, dhi_req);

if ~(isfinite(dlo) && isfinite(dhi) && dhi > dlo)
    error("Requested scan window [%.6g, %.6g] has no overlap with data range [%.6g, %.6g].", ...
        dlo_req, dhi_req, dmin_all, dmax_all);
end

if mode == "step"
    step = double(hyst.step);
    if ~isfinite(step) || step <= 0
        error("hyst.step must be > 0 when grid_mode='step'");
    end
    dop_base = (dlo:step:dhi).';
    if isempty(dop_base) || abs(dop_base(end)-dhi) > 1e-12
        dop_base(end+1,1) = dhi;
    end
else
    N = double(hyst.N);
    if ~isfinite(N) || N < 2
        error("hyst.N must be >= 2 when grid_mode='N'");
    end
    dop_base = linspace(dlo, dhi, round(N)).';
end

if direction == "ascend"
    dop_fwd = dop_base;
else
    dop_fwd = flipud(dop_base);
end

dop_bwd = flipud(dop_fwd);
end

function [Yf_plot, Yb_plot, gamma_used] = apply_lorentz_on_doping_(xf, Yf, xb, Yb, smooth)
% plotting-only smoothing
Yf_plot = Yf;
Yb_plot = Yb;
gamma_used = NaN;

if ~isstruct(smooth), return; end
if ~isfield(smooth, 'use_lorentz') || ~smooth.use_lorentz, return; end
if ~isfield(smooth, 'gamma') || ~isfinite(smooth.gamma) || smooth.gamma <= 0, return; end

gamma_used = smooth.gamma;
Yf_plot = lorentz_smooth_1d_(xf, Yf, gamma_used);
Yb_plot = lorentz_smooth_1d_(xb, Yb, gamma_used);
end

function ys = lorentz_smooth_1d_(x, y, gamma)
x = x(:);
y = y(:);
ys = nan(size(y));

for i = 1:numel(x)
    w = gamma^2 ./ ((x - x(i)).^2 + gamma^2);
    ok = isfinite(w) & isfinite(y);
    if any(ok)
        ys(i) = sum(w(ok).*y(ok)) / sum(w(ok));
    end
end
end