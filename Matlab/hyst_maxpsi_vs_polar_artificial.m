function out = hyst_maxpsi_vs_polar_artificial()
% hyst_maxpsi_vs_polar_artificial
% Use polar=0 baseline chi0(mu) to build artificial polar-field chi_p(mu) by shift-average:
%   chi_p(mu) = 0.5*(chi0(mu+dmu)+chi0(mu-dmu)), dmu = 0.5*polar_meV*(meV->mu_unit)
% Run hysteresis (forward + backward) on valid mu ticks in par.io.mu_range,
% then compute max|psi| within that mu-range vs polar.
%
% NOTE: Minimal "sanity checks": only what prevents crashes.
%       Added progress + debug prints + optional waitbar.

% ---------------- user settings ----------------
T_target = 4;          % K
T_tol    = 1e-6;       % K
iq_pick  = -133;
jq_pick  = -133;

polar_list_meV    = linspace(0, 3.6, 19);  % meV
mu_unit_in_header = "eV";                  % "eV" or "meV"

use_branch   = "forward";  % "forward" | "backward" | "max_of_both"
use_dop_range = false;     % optionally also restrict by par.io.dop_range

% progress/debug
dbg = struct();
dbg.enable          = true;
dbg.use_waitbar     = true;
dbg.print_every_mu  = 10;     % print every N mu ticks
dbg.print_every_file= 200;    % baseline scan print every N files
dbg.slow_sec        = 0.75;   % warn if a step slower than this
% -----------------------------------------------

par  = make_realchi_params();
coef = make_realchi_coeff(par);

% range in header mu unit
mu_lo = min(par.io.mu_range);
mu_hi = max(par.io.mu_range);

% choose baseline root (polar=0 dataset)
default_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_200\D0.084";
if ~isfolder(default_root), default_root = string(pwd); end
root = uigetdir(default_root, "Select BASELINE (polar=0) root folder (recursive chi*.txt)");
if isequal(root,0), out = struct(); return; end
root = string(root);

% meV -> header mu shift unit
switch lower(string(mu_unit_in_header))
    case "ev"
        polar_to_mu = 1e-3;  % meV -> eV
    case "mev"
        polar_to_mu = 1.0;   % meV -> meV
    otherwise
        error('mu_unit_in_header must be "eV" or "meV".');
end

% optional waitbar
wb = [];
if dbg.enable && dbg.use_waitbar
    wb = waitbar(0, "Starting...", "Name", "hyst_maxpsi_vs_polar");
end
cleanupWb = onCleanup(@() safe_close_waitbar_(wb)); %#ok<NASGU>

% 1) baseline: (mu, dop, chi)
if dbg.enable
    fprintf("=== Collect baseline (polar=0) ===\n");
end
B = collect_baseline_curve_(root, T_target, T_tol, iq_pick, jq_pick, dbg, wb);

mu0  = B.mu(:);
dop0 = B.dop(:);
chi0 = B.chi(:);

% sort by mu
[mu0, ord] = sort(mu0);
dop0 = dop0(ord);
chi0 = chi0(ord);

% interpolants
Fchi0 = griddedInterpolant(mu0, chi0, "linear", "none");
Fdop0 = griddedInterpolant(mu0, dop0, "linear", "none");

% pack scan/opt from par.min
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

P = unique(double(polar_list_meV(:)));
P = P(isfinite(P));
P = sort(P);

max_abspsi = nan(size(P));
mu_at_max  = nan(size(P));
dop_at_max = nan(size(P));

if dbg.enable
    fprintf("=== Hysteresis + max|psi| vs polar ===\n");
    fprintf("mu-range=[%.6g, %.6g] (%s)\n", mu_lo, mu_hi, mu_unit_in_header);
    fprintf("branch=%s, q=(%d,%d), T=%.6gK\n", use_branch, iq_pick, jq_pick, T_target);
end

for ip = 1:numel(P)
    pmev = P(ip);
    dmu  = 0.5 * pmev * polar_to_mu;

    if dbg.enable
        fprintf("\n--- polar %d/%d : %.6g meV (dmu=%.6g %s) ---\n", ...
            ip, numel(P), pmev, dmu, mu_unit_in_header);
    end
    if ~isempty(wb) && isvalid(wb)
        waitbar((ip-1)/numel(P), wb, sprintf("polar %d/%d : %.3g meV", ip, numel(P), pmev));
        drawnow;
    end

    % shift-average chi_p on baseline ticks, then keep finite
    chi_p = 0.5 * (Fchi0(mu0 + dmu) + Fchi0(mu0 - dmu));
    ok = isfinite(chi_p);
    if nnz(ok) < 5
        if dbg.enable, fprintf("skip: too few finite points after shift-average (%d)\n", nnz(ok)); end
        continue;
    end

    mu_scan  = mu0(ok);
    chi_scan = chi_p(ok);
    dop_scan = Fdop0(mu_scan);

    % collapse duplicates
    [mu_scan, chi_scan, dop_scan] = collapse_same_mu_(mu_scan, chi_scan, dop_scan);
    if numel(mu_scan) < 5
        if dbg.enable, fprintf("skip: too few points after collapse (%d)\n", numel(mu_scan)); end
        continue;
    end

    % restrict mu-range
    inMu = (mu_scan >= mu_lo) & (mu_scan <= mu_hi);
    if nnz(inMu) < 5
        if dbg.enable
            fprintf("skip: mu in range has only %d points\n", nnz(inMu));
        end
        continue;
    end
    mu_scan  = mu_scan(inMu);
    chi_scan = chi_scan(inMu);
    dop_scan = dop_scan(inMu);

    % optional dop-range
    if use_dop_range && isfield(par,'io') && isfield(par.io,'dop_range')
        dlo = min(par.io.dop_range);
        dhi = max(par.io.dop_range);
        inD = (dop_scan >= dlo) & (dop_scan <= dhi);
        if nnz(inD) < 5
            if dbg.enable, fprintf("skip: dop in range has only %d points\n", nnz(inD)); end
            continue;
        end
        mu_scan  = mu_scan(inD);
        chi_scan = chi_scan(inD);
        dop_scan = dop_scan(inD);
    end

    % sort increasing mu
    [mu_f, ord2] = sort(mu_scan(:));
    chi_scan = chi_scan(ord2);
    dop_scan = dop_scan(ord2);

    mu_b = flipud(mu_f);

    Fchi = griddedInterpolant(mu_f, chi_scan(:), "linear", "none");
    Fdop = griddedInterpolant(mu_f, dop_scan(:), "linear", "none");

    % forward
    psi_f = nan(numel(mu_f),2);
    prev  = [par.hyst.psi1_0, par.hyst.psi2_0];
    tPolar = tic;

    for i = 1:numel(mu_f)
        mu  = mu_f(i);
        dop = Fdop(mu);
        chi_used = Fchi(mu);

        C = coef.eval(T_target, dop, chi_used);
        F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);

        tStep = tic;
        mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);
        dt = toc(tStep);

        if dbg.enable && dt > dbg.slow_sec
            fprintf("  [slow] FWD i=%d/%d mu=%.6g dop=%.6g mins=%.3fs\n", ...
                i, numel(mu_f), mu, dop, dt);
        end

        if isempty(mins)
            [p1,p2] = minimize_free_energy(C.a1,C.b1,C.a2,C.b2,C.c2,par.lambda, prev(1),prev(2), false);
            pick = [p1,p2];
        else
            pick = pick_closest_min_(mins, prev);
        end

        psi_f(i,:) = pick;
        prev = pick;

        if dbg.enable && (mod(i, dbg.print_every_mu)==0 || i==1 || i==numel(mu_f))
            elapsed = toc(tPolar);
            frac = i/numel(mu_f);
            eta = elapsed*(1/frac - 1);
            fprintf("  [FWD] %4d/%4d (%.1f%%) elapsed=%.1fs ETA=%.1fs\n", ...
                i, numel(mu_f), 100*frac, elapsed, eta);
            if ~isempty(wb) && isvalid(wb)
                waitbar((ip-1 + 0.35*frac)/numel(P), wb, ...
                    sprintf("polar %.3g meV FWD %d/%d", pmev, i, numel(mu_f)));
            end
            drawnow;
        end
    end

    % backward
    psi_b = nan(numel(mu_f),2);
    prev  = psi_f(end,:);
    for i = 1:numel(mu_b)
        mu  = mu_b(i);
        dop = Fdop(mu);
        chi_used = Fchi(mu);

        C = coef.eval(T_target, dop, chi_used);
        F = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, par.lambda);

        tStep = tic;
        mins = find_local_minima_2d_denseSeeds_(F, par.psi1_lim, par.psi2_lim, scan, opt);
        dt = toc(tStep);

        if dbg.enable && dt > dbg.slow_sec
            fprintf("  [slow] BWD i=%d/%d mu=%.6g dop=%.6g mins=%.3fs\n", ...
                i, numel(mu_b), mu, dop, dt);
        end

        if isempty(mins)
            [p1,p2] = minimize_free_energy(C.a1,C.b1,C.a2,C.b2,C.c2,par.lambda, prev(1),prev(2), false);
            pick = [p1,p2];
        else
            pick = pick_closest_min_(mins, prev);
        end

        psi_b(i,:) = pick;
        prev = pick;

        if dbg.enable && (mod(i, dbg.print_every_mu)==0 || i==1 || i==numel(mu_b))
            elapsed = toc(tPolar);
            frac = i/numel(mu_b);
            eta = elapsed*(1/frac - 1);
            fprintf("  [BWD] %4d/%4d (%.1f%%) elapsed=%.1fs ETA=%.1fs\n", ...
                i, numel(mu_b), 100*frac, elapsed, eta);
            if ~isempty(wb) && isvalid(wb)
                waitbar((ip-1 + 0.70 + 0.25*frac)/numel(P), wb, ...
                    sprintf("polar %.3g meV BWD %d/%d", pmev, i, numel(mu_b)));
            end
            drawnow;
        end
    end

    psi_b_aligned = flipud(psi_b);

    % max|psi|
    abspsi_f = abs(psi_f(:,1));
    abspsi_b = abs(psi_b_aligned(:,1));

    switch lower(string(use_branch))
        case "forward"
            [v, ix] = max(abspsi_f);
        case "backward"
            [v, ix] = max(abspsi_b);
        case "max_of_both"
            [v1, ix1] = max(abspsi_f);
            [v2, ix2] = max(abspsi_b);
            if v2 > v1
                v = v2; ix = ix2;
            else
                v = v1; ix = ix1;
            end
        otherwise
            error('use_branch must be "forward" | "backward" | "max_of_both".');
    end

    max_abspsi(ip) = v;
    mu_at_max(ip)  = mu_f(ix);
    dop_at_max(ip) = Fdop(mu_f(ix));

    if dbg.enable
        fprintf("  ==> max|psi|=%.6g at mu=%.6g dop=%.6g\n", ...
            max_abspsi(ip), mu_at_max(ip), dop_at_max(ip));
    end
end

% plot
fig = figure("Color","w","Units","pixels","Position",[120 120 900 420], ...
    "Name","max|psi| vs polar (artificial shift-average)");
ax = axes(fig);
plot(ax, P, max_abspsi, "o-","LineWidth",1.8,"MarkerSize",6);
grid(ax,"on"); box(ax,"on");
xlabel(ax,"polar\_mu (meV)");
ylabel(ax,"max |psi| within mu-range");
title(ax, sprintf("T=%.6gK q=(%d,%d) mu-range=[%.6g,%.6g] (%s) branch=%s", ...
    T_target, iq_pick, jq_pick, mu_lo, mu_hi, mu_unit_in_header, use_branch), ...
    "FontWeight","normal");

out = struct();
out.root = root;
out.T_target = T_target;
out.iq_pick = iq_pick;
out.jq_pick = jq_pick;
out.mu_unit_in_header = mu_unit_in_header;
out.mu_range = [mu_lo, mu_hi];
out.use_branch = use_branch;
out.polar_meV = P;
out.max_abspsi = max_abspsi;
out.mu_at_max  = mu_at_max;
out.dop_at_max = dop_at_max;

if ~isempty(wb) && isvalid(wb)
    waitbar(1, wb, "Done.");
    drawnow;
end
end

% ============================================================
% helpers
% ============================================================

function B = collect_baseline_curve_(root, T_target, T_tol, iq_pick, jq_pick, dbg, wb)
if nargin < 6 || isempty(dbg), dbg = struct('enable',false); end
if nargin < 7, wb = []; end

L = dir(fullfile(root, "**", "chi*.txt"));
L = L(~[L.isdir]);

Tlo = T_target - T_tol;
Thi = T_target + T_tol;

mu  = nan(numel(L),1);
dop = nan(numel(L),1);
chi = nan(numel(L),1);

n=0; dropT=0; dropH=0; dropTbl=0; dropQ=0;

for k = 1:numel(L)
    if dbg.enable && (mod(k, dbg.print_every_file)==0 || k==1 || k==numel(L))
        fprintf("[baseline scan] %d/%d (kept=%d)\n", k, numel(L), n);
        if ~isempty(wb) && isvalid(wb)
            waitbar(min(0.05 + 0.20*(k/numel(L)), 0.25), wb, ...
                sprintf("Scanning baseline %d/%d (kept=%d)", k, numel(L), n));
            drawnow;
        end
    end

    fpath = string(fullfile(L(k).folder, L(k).name));

    H = parse_header_T_mu_dop_(fpath);

    if ~isfinite(H.T), dropH=dropH+1; continue; end
    if ~(H.T>=Tlo && H.T<=Thi), dropT=dropT+1; continue; end

    if ~isfinite(H.mu) || ~isfinite(H.doping)
        [mu_p, dop_p] = parse_mu_dop_from_path_(string(L(k).folder));
        if ~isfinite(H.mu), H.mu = mu_p; end
        if ~isfinite(H.doping), H.doping = dop_p; end
    end
    if ~isfinite(H.mu) || ~isfinite(H.doping), dropH=dropH+1; continue; end

    try
        M = readmatrix(fpath, "FileType","text", "CommentStyle","#");
    catch
        dropTbl=dropTbl+1; continue;
    end
    if isempty(M), dropTbl=dropTbl+1; continue; end

    Cc = detect_cols_chi_(M);
    id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
    if isempty(id), dropQ=dropQ+1; continue; end

    n=n+1;
    mu(n)  = H.mu;
    dop(n) = H.doping;
    chi(n) = M(id, Cc.Re);
end

mu  = mu(1:n);
dop = dop(1:n);
chi = chi(1:n);

[mu_u,~,ic] = unique(double(mu(:)));
chi_m = accumarray(ic, double(chi(:)), [], @mean);
dop_m = accumarray(ic, double(dop(:)), [], @mean);

B = struct();
B.mu  = mu_u;
B.dop = dop_m;
B.chi = chi_m;

if dbg.enable
    fprintf("[baseline scan done] files=%d kept=%d dropT=%d dropH=%d dropTbl=%d dropQ=%d\n", ...
        numel(L), n, dropT, dropH, dropTbl, dropQ);
end
end

function H = parse_header_T_mu_dop_(fpath)
H = struct('T',NaN,'mu',NaN,'doping',NaN);

fid = fopen(fpath,'r');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    s0 = strtrim(string(ln));
    if ~startsWith(s0,"#"), break; end
    s = strtrim(erase(s0,"#"));

    if ~isfinite(H.T)
        tok = regexp(s, "(?:^|\s)T(?:_K)?\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.T = str2double(tok{1}); end
    end
    if ~isfinite(H.doping)
        tok = regexp(s, "(?:^|\s)(?:doping|dop)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.doping = str2double(tok{1}); end
    end
    if ~isfinite(H.mu)
        tok = regexp(s, "(?:^|\s)(?:mu|mu_eV|EF|EF_eV|E_F)\s*=\s*" + num, "tokens","once");
        if ~isempty(tok), H.mu = str2double(tok{1}); end
    end
end
end

function [mu_val, dop_val] = parse_mu_dop_from_path_(folder)
mu_val = NaN; dop_val = NaN;
p = replace(string(folder),"\","/");

tok = regexp(p, "/mu([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)(?:/|$)", "tokens","once","ignorecase");
if ~isempty(tok), mu_val = str2double(tok{1}); end

tok = regexp(p, "/doping([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)(?:/|$)", "tokens","once","ignorecase");
if ~isempty(tok), dop_val = str2double(tok{1}); end
end

function C = detect_cols_chi_(M)
% supports either:
%   >=8 cols: idx iq jq qx qy Re Im ...
%   else:     iq jq qx qy Re Im ...
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
[mu_u,~,ic] = unique(mu);
chi_m = accumarray(ic, chi, [], @mean);
dop_m = accumarray(ic, dop, [], @mean);
[mu2,ord] = sort(mu_u);
chi2 = chi_m(ord);
dop2 = dop_m(ord);
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
    lin = sub2ind(size(P1), rr, cc);
    seeds = [P1(lin), P2(lin)];
    E = Fgrid(lin);
    [~,ix] = sort(E,'ascend');
    ix = ix(1:min(scan.keep_topK, numel(ix)));
    seeds = seeds(ix,:);
end

if opt.force_origin_and_nearby
    d = opt.nearby_delta;
    seeds = [seeds; 0,0; d,0; -d,0; 0,d; 0,-d];
end

obj = @(x) safe_F_(F, x(1), x(2), psi1_lim, psi2_lim);
fopt = optimset('Display','off', 'MaxIter', opt.max_iter, 'TolFun', opt.tol_fun, 'TolX', opt.tol_x);

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
    Hh = hessian_fd_(@(u,v) F(u,v), x(1), x(2), opt.fd_h);
    ev = eig((Hh+Hh')/2);
    if all(real(ev) > opt.min_eig_eps)
        keep(i) = true;
    end
end
mins = uniq(keep,:);
end

function val = safe_F_(F, psi1, psi2, lim1, lim2)
if psi1 < lim1(1) || psi1 > lim1(2) || psi2 < lim2(1) || psi2 > lim2(2)
    val = 1e30; return;
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
n1 = A(1:R-2, 2:C-1); n2 = A(3:R, 2:C-1);
n3 = A(2:R-1, 1:C-2); n4 = A(2:R-1, 3:C);
n5 = A(1:R-2, 1:C-2); n6 = A(1:R-2, 3:C);
n7 = A(3:R, 1:C-2);   n8 = A(3:R, 3:C);
m = center <= n1 & center <= n2 & center <= n3 & center <= n4 & ...
    center <= n5 & center <= n6 & center <= n7 & center <= n8;
mask(2:R-1, 2:C-1) = m;
end

function B = gaussian_smooth_2d_(A, sigma)
if sigma <= 0, B = A; return; end
rad = max(1, ceil(3*sigma));
x = (-rad:rad);
g = exp(-(x.^2)/(2*sigma^2)); g = g / sum(g);
B = conv2(A, g, 'same');
B = conv2(B, g', 'same');
end

function safe_close_waitbar_(wb)
try
    if ~isempty(wb) && isvalid(wb)
        close(wb);
    end
catch
end
end