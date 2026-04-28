function OUT = hyst_maxpsi_vs_polar_dense()
% hyst_maxpsi_vs_polar_dense
% ------------------------------------------------------------
% Read chi data ONLY ONCE, then scan custom artificial polar list.
% For each polar:
%   1) build chi_art(mu) = 0.5*(chi(mu+dmu)+chi(mu-dmu))
%   2) map to doping
%   3) run forward/backward hysteresis continuation
%      using dense-seed local minima search
%   4) extract max |psi1|
%
% Based on:
%   - hysteresis_realchi_T_SinglePolar_LorentzBroadening
%   - hyst_maxpsi_vs_polar_cached
% ------------------------------------------------------------

% ---------------- params / coeff ----------------
par  = make_realchi_params(true);
coef = make_realchi_coeff(par);

% ---------------- choose folder -----------------
default_root = "E:/rg_master/data/";
if ~isfolder(default_root)
    default_root = string(pwd);
end
root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
if isequal(root,0)
    error('User cancelled.');
end
root = string(root);

fprintf("Root: %s\n", root);

% ---------------- custom polar list -------------
polarList_meV = 0:0.2:8;   % <<< 改这里

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

if numel(mu0) < 5
    error("Too few points at this T. Check T_tol / header parse / (iq,jq).");
end

fprintf("[baseline] points=%d | mu=[%.6g, %.6g] | dop=[%.6g, %.6g]\n", ...
    numel(mu0), min(mu0), max(mu0), min(dop0), max(dop0));

% sort by mu & collapse duplicates
[mu0, chi0, dop0] = collapse_same_mu_(mu0, chi0, dop0);

% interpolants on mu
Fchi0 = griddedInterpolant(mu0, chi0, "linear", "none");
Fdop0 = griddedInterpolant(mu0, dop0, "linear", "none");

% common doping path once
[dop_fwd, dop_bwd] = build_hyst_doping_grid_(dop0, par.hyst);
Nscan = numel(dop_fwd);

fprintf("[scan window] forward: [%.6g, %.6g], backward: [%.6g, %.6g], N=%d\n", ...
    dop_fwd(1), dop_fwd(end), dop_bwd(1), dop_bwd(end), Nscan);

% ---------------- pack scan/opt from par.min ----
scan = struct();
scan.Npsi1 = par.min.Npsi1;
scan.Npsi2 = par.min.Npsi2;
scan.use_gaussian_smooth = par.min.use_gaussian_smooth;
scan.smooth_sigma = par.min.smooth_sigma;
scan.keep_topK = par.min.keep_topK;

opt = struct();
opt.seed_jitter  = par.min.seed_jitter;
opt.seed_repeats = par.min.seed_repeats;
opt.max_iter     = par.min.max_iter;
opt.tol_fun      = par.min.tol_fun;
opt.tol_x        = par.min.tol_x;
opt.cluster_tol  = par.min.cluster_tol;
opt.fd_h         = par.min.fd_h;
opt.min_eig_eps  = par.min.min_eig_eps;
opt.force_origin_and_nearby = par.min.force_origin_and_nearby;
opt.nearby_delta = par.min.nearby_delta;

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
    % forward continuation (dense minima version)
    % -----------------------------------------------------
    psi_fwd = nan(Nscan,2);
    a1_fwd  = nan(Nscan,1);
    chi_fwd = nan(Nscan,1);

    prev = [par.hyst.psi1_0, par.hyst.psi2_0];

    for i = 1:Nscan
        dop = dop_fwd(i);
        chi_now = chi_of_dop(dop);
        C = coef.eval(T_target, dop, chi_now);

        F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
        mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);

        if isempty(mins)
            [p1,p2,~] = minimize_free_energy( ...
                C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda, ...
                prev(1), prev(2), false);
            pick = [p1,p2];
        else
            pick = pick_closest_min_(mins, prev);
        end

        psi_fwd(i,:) = pick;
        prev = pick;

        a1_fwd(i)  = C.a1;
        chi_fwd(i) = chi_now;
    end

    % -----------------------------------------------------
    % backward continuation (dense minima version)
    % -----------------------------------------------------
    psi_bwd = nan(Nscan,2);
    a1_bwd  = nan(Nscan,1);
    chi_bwd = nan(Nscan,1);

    prev = psi_fwd(end,:);

    for i = 1:Nscan
        dop = dop_bwd(i);
        chi_now = chi_of_dop(dop);
        C = coef.eval(T_target, dop, chi_now);

        F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);
        mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);

        if isempty(mins)
            [p1,p2,~] = minimize_free_energy( ...
                C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda, ...
                prev(1), prev(2), false);
            pick = [p1,p2];
        else
            pick = pick_closest_min_(mins, prev);
        end

        psi_bwd(i,:) = pick;
        prev = pick;

        a1_bwd(i)  = C.a1;
        chi_bwd(i) = chi_now;
    end

    % -----------------------------------------------------
    % optional plotting-only Lorentz smoothing on psi1
    % -----------------------------------------------------
    psi_f_abs = abs(psi_fwd);
    psi_b_abs = abs(psi_bwd);

    [psi_f_abs_plot, psi_b_abs_plot, gamma_used] = apply_lorentz_on_doping_( ...
        dop_fwd, psi_f_abs, dop_bwd, psi_b_abs, par.smooth);

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
    H.polar_meV   = pmev;
    H.dmu_eV      = dmu;
    H.mu_used     = mu_use;
    H.dop_used    = dop_use;
    H.chi_used    = chi_use;
    H.dop_fwd     = dop_fwd;
    H.dop_bwd     = dop_bwd;
    H.psi_fwd     = psi_fwd;
    H.psi_bwd     = psi_bwd;
    H.psi1_f_plot = psi_f_abs_plot(:,1);
    H.psi1_b_plot = psi_b_abs_plot(:,1);
    H.gamma_used  = gamma_used;
    H.a1_fwd      = a1_fwd;
    H.a1_bwd      = a1_bwd;
    H.chi_fwd     = chi_fwd;
    H.chi_bwd     = chi_bwd;
    H.maxPsi_fwd  = maxPsi_fwd(ip);
    H.maxPsi_bwd  = maxPsi_bwd(ip);
    H.maxPsi_all  = maxPsi_all(ip);

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
    "maxPsi_vs_polar_dense_T%.4f_iq%d_jq%d.png", ...
    T_target, iq_pick, jq_pick));

out_mat = fullfile(out_dir, sprintf( ...
    "maxPsi_vs_polar_dense_T%.4f_iq%d_jq%d.mat", ...
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
if isempty(L), error("No chi*.txt under: %s", root); end

Tlo = T_target - T_tol;
Thi = T_target + T_tol;

mu  = nan(numel(L),1);
dop = nan(numel(L),1);
chi = nan(numel(L),1);

n=0; dropT=0; dropH=0; dropTbl=0; dropQ=0;

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
    if ~isfinite(T_use), dropH=dropH+1; continue; end
    if ~(T_use>=Tlo && T_use<=Thi), dropT=dropT+1; continue; end

    mu_path = parse_mu_from_path_(folder);
    mu_use = mu_path;
    if ~isfinite(mu_use), mu_use = H.mu; end

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
    id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
    if isempty(id), dropQ=dropQ+1; continue; end

    n=n+1;
    mu(n)  = mu_use;
    dop(n) = H.dop;
    chi(n) = M(id, Cc.Re);
end

mu=mu(1:n); dop=dop(1:n); chi=chi(1:n);
ok = isfinite(mu) & isfinite(dop) & isfinite(chi);
mu=mu(ok); dop=dop(ok); chi=chi(ok);

fprintf("[collect] files=%d kept=%d | dropT=%d dropH=%d dropTbl=%d dropQ=%d\n", ...
    numel(L), numel(mu), dropT, dropH, dropTbl, dropQ);

B = struct('mu',mu,'dop',dop,'chi',chi);
end

function H = parse_header_mu_dop_T_(fpath)
H = struct('T',NaN,'mu',NaN,'dop',NaN);
fid = fopen(fpath,'r'); if fid<0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num  = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

for t = 1:3000
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    s0 = strtrim(string(ln));
    if ~startsWith(s0,"#"), break; end
    s = strtrim(erase(s0,"#"));

    if ~isfinite(H.T)
        tok = regexp(s, "(?:^|\s)T(?:_K)?\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.T = str2double(tok{1}); end
    end
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

function T = parse_T_from_path_(folder)
folder = replace(string(folder),"\","/");
num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

tok = regexp(folder, "/T[_=]?\s*" + num, "tokens","once","ignorecase");
if ~isempty(tok)
    T = str2double(tok{1});
    return;
end

tok = regexp(folder, "/temp[_=]?\s*" + num, "tokens","once","ignorecase");
if ~isempty(tok)
    T = str2double(tok{1});
else
    T = NaN;
end
end

function mu = parse_mu_from_path_(folder)
folder = replace(string(folder),"\","/");
num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";
tok = regexp(folder, "/mu\s*=?\s*" + num + "(?:/|$)", "tokens","once","ignorecase");
if ~isempty(tok), mu = str2double(tok{1}); else, mu = NaN; end
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
mu=mu(ok); chi=chi(ok); dop=dop(ok);

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
dop=dop(ok); chi=chi(ok); mu=mu(ok);

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
if isfield(hyst, "grid_mode") && ~isempty(hyst.grid_mode)
    mode = string(hyst.grid_mode);
end

direction = "ascend";
if isfield(hyst, "forward_direction") && ~isempty(hyst.forward_direction)
    direction = lower(string(hyst.forward_direction));
end

if ~(direction == "ascend" || direction == "descend")
    error("hyst.forward_direction must be 'ascend' or 'descend'");
end

scan_start = NaN;
scan_end   = NaN;

if isfield(hyst, "scan_start") && ~isempty(hyst.scan_start)
    scan_start = double(hyst.scan_start);
end
if isfield(hyst, "scan_end") && ~isempty(hyst.scan_end)
    scan_end = double(hyst.scan_end);
end

if ~isfinite(scan_start)
    scan_start = dmin_all;
end
if ~isfinite(scan_end)
    scan_end = dmax_all;
end

dlo_req = min(scan_start, scan_end);
dhi_req = max(scan_start, scan_end);

dlo = max(dmin_all, dlo_req);
dhi = min(dmax_all, dhi_req);

if ~(isfinite(dlo) && isfinite(dhi) && dhi > dlo)
    error("Requested scan window [%.6g, %.6g] has no overlap with data range [%.6g, %.6g].", ...
        dlo_req, dhi_req, dmin_all, dmax_all);
end

fprintf("[build grid] data range = [%.6g, %.6g], requested = [%.6g, %.6g], used = [%.6g, %.6g]\n", ...
    dmin_all, dmax_all, dlo_req, dhi_req, dlo, dhi);

if mode == "step"
    step = double(hyst.step);
    if ~isfinite(step) || step <= 0
        error("hyst.step must be > 0 when grid_mode='step'");
    end

    dop_base = (dlo:step:dhi).';
    if isempty(dop_base) || abs(dop_base(end) - dhi) > 1e-12
        dop_base(end+1,1) = dhi; %#ok<AGROW>
    end

elseif mode == "N"
    N = double(hyst.N);
    if ~isfinite(N) || N < 2
        error("hyst.N must be >= 2 when grid_mode='N'");
    end

    dop_base = linspace(dlo, dhi, round(N)).';

else
    error("Unsupported hyst.grid_mode. Use 'N' or 'step'.");
end

if direction == "ascend"
    dop_fwd = dop_base;
else
    dop_fwd = flipud(dop_base);
end

dop_bwd = flipud(dop_fwd);
end

function [Yf_plot, Yb_plot, gamma_used] = apply_lorentz_on_doping_(xf, Yf, xb, Yb, smooth)
Yf_plot = Yf;
Yb_plot = Yb;
gamma_used = NaN;

if ~isfield(smooth,'use_lorentz') || ~smooth.use_lorentz
    return;
end

dx = median(abs(diff(xf)));
if ~isfinite(dx) || dx <= 0
    return;
end

gfac = 0;
if isfield(smooth,'lorentz_gamma_factor')
    gfac = smooth.lorentz_gamma_factor;
end

if ~isfinite(gfac) || gfac < 0
    return;
end

if gfac == 0
    gamma = 1.5 * dx;
else
    gamma = gfac * dx;
end

if ~isfinite(gamma) || gamma <= 0
    return;
end
gamma_used = gamma;

Yf_plot = zeros(size(Yf));
Yb_plot = zeros(size(Yb));
for c = 1:size(Yf,2)
    Yf_plot(:,c) = lorentz_smooth_1d_(xf, Yf(:,c), gamma);
end
for c = 1:size(Yb,2)
    Yb_plot(:,c) = lorentz_smooth_1d_(xb, Yb(:,c), gamma);
end
end

function y_s = lorentz_smooth_1d_(x, y, gamma)
x = x(:); y = y(:);
N = numel(x);
y_s = zeros(N,1);

for i = 1:N
    dx = x(i) - x;
    w  = (gamma/pi) ./ (dx.^2 + gamma^2);
    sw = sum(w);
    if sw > 0
        y_s(i) = sum(w .* y) / sw;
    else
        y_s(i) = y(i);
    end
end
end

function pick = pick_closest_min_(mins, prev)
d = sqrt(sum((mins - prev).^2, 2));
[~,ix] = min(d);
pick = mins(ix,:);
end

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

if isfield(opt,'force_origin_and_nearby') && opt.force_origin_and_nearby
    d = opt.nearby_delta;
    seeds = [seeds;
             0,0;
             d,0; -d,0; 0,d; 0,-d];
end

obj = @(x) safe_F_(F, x(1), x(2), psi1_lim, psi2_lim);

fopt = optimset('Display','off', ...
    'MaxIter', opt.max_iter, ...
    'TolFun',  opt.tol_fun, ...
    'TolX',    opt.tol_x);

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
    Hh = hessian_fd_( @(u,v) F(u,v), x(1), x(2), opt.fd_h );
    ev = eig((Hh+Hh')/2);
    if all(real(ev) > opt.min_eig_eps)
        keep(i) = true;
    end
end

mins = uniq(keep,:);
end

function val = safe_F_(F, psi1, psi2, lim1, lim2)
if psi1 < lim1(1) || psi1 > lim1(2) || psi2 < lim2(1) || psi2 > lim2(2)
    val = 1e30;
    return;
end
val = F(psi1, psi2);
if ~isfinite(val)
    val = 1e30;
end
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

function Hh = hessian_fd_(f, x, y, h)
f00 = f(x,y);
fxx = (f(x+h,y) - 2*f00 + f(x-h,y))/h^2;
fyy = (f(x,y+h) - 2*f00 + f(x,y-h))/h^2;
fxy = (f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)) / (4*h^2);
Hh = [fxx, fxy; fxy, fyy];
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
if sigma <= 0
    B = A;
    return;
end
rad = max(1, ceil(3*sigma));
x = (-rad:rad);
g = exp(-(x.^2)/(2*sigma^2));
g = g / sum(g);
B = conv2(A, g, 'same');
B = conv2(B, g', 'same');
end