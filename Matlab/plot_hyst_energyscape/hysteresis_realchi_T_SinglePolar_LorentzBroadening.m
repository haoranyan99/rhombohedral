 function H = hysteresis_realchi_T_SinglePolar_LorentzBroadening()
% Hysteresis at ONE temperature using chi*.txt (recursive).
% IMPORTANT: chi panel + coefficients a are computed from SHIFT-AVERAGED chi:
%   chi_art(mu) = 0.5*(chi(mu+dmu) + chi(mu-dmu)),
%   dmu = 0.5 * artificial_polar * 1e-3 (eV).
%
% Needs existing functions:
%   - make_realchi_params
%   - make_realchi_coeff
%   - free_energy
%   - minimize_free_energy

% ---------------- params / coeff ----------------
par  = make_realchi_params(true);
coef = make_realchi_coeff(par);


% ---------------- choose folder -----------------
% default_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_200\D0.084\polar_meV0.000";
default_root = "E:\rg_master\data";
if ~isfolder(default_root)
    default_root = string(pwd);
end
root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
if isequal(root, 0), error('User cancelled.'); end
root = string(root);

fprintf("Root: %s\n", root);

% ---------------- q point (ABS) -----------------
iq_pick = round(par.iq_pick);
jq_pick = round(par.jq_pick);
fprintf("Pick (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);

% ---------------- target T -----------------
T_target = par.T_target;
T_tol    = par.T_tol;
fprintf("T_target=%.6g K, T_tol=%.3g K\n", T_target, T_tol);

% ---------------- artificial polar shift ----------------
pmev = double(par.artificial_polar);
dmu  = 0.5 * pmev * 1e-3;   % meV -> eV, then /2
fprintf("artificial_polar=%.6g meV  => dmu=%.6g eV\n", pmev, dmu);

% =========================================================
% 1) collect baseline curve at this T: (mu_i, doping_i, chi_i)
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

% sort by mu & collapse duplicates
[mu0, chi0, dop0] = collapse_same_mu_(mu0, chi0, dop0);

% build interpolants on mu
Fchi0 = griddedInterpolant(mu0, chi0, "linear", "none");
Fdop0 = griddedInterpolant(mu0, dop0, "linear", "none");

% =========================================================
% 2) build chi_art(mu) on mu ticks, then map to doping axis
% =========================================================
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

% drop non-finite dop too
ok2 = isfinite(mu_use) & isfinite(chi_use) & isfinite(dop_use);
mu_use  = mu_use(ok2);
chi_use = chi_use(ok2);
dop_use = dop_use(ok2);

if numel(mu_use) < 5
    error("After shift-average, too few finite points. (dmu too large or mu range too narrow)");
end

% sort by doping (because hyst uses doping axis)
[dop_use, ord] = sort(dop_use);
chi_use = chi_use(ord);
mu_use  = mu_use(ord);

% collapse duplicate doping (mean)
[dop_use, chi_use, mu_use] = collapse_same_dop_(dop_use, chi_use, mu_use);

% interpolation chi_art(dop)
chi_of_dop = @(d) interp1(dop_use, chi_use, d, 'linear', 'extrap');

% =========================================================
% 3) build hyst doping grids
% =========================================================
[dop_fwd, dop_bwd] = build_hyst_doping_grid_(dop_use, par.hyst);
N = numel(dop_fwd);

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
% 4) forward continuation (USES chi_art!)
% =========================================================
psi_f = nan(N,2);
a1_f  = nan(N,1);
chi_f = nan(N,1);

prev = [par.hyst.psi1_0, par.hyst.psi2_0];

for i = 1:N
    dop = dop_fwd(i);
    chi_used = chi_of_dop(dop);              % <-- chi_art(dop)
    C = coef.eval(T_target, dop, chi_used);  % <-- coefficients from chi_art

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

% =========================================================
% 5) backward continuation (USES chi_art!)
% =========================================================
psi_b = nan(N,2);
a1_b  = nan(N,1);
chi_b = nan(N,1);

prev = psi_f(end,:);

for i = 1:N
    dop = dop_bwd(i);
    chi_used = chi_of_dop(dop);
    C = coef.eval(T_target, dop, chi_used);

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

% =========================================================
% 6) Lorentz broadening for plotting (only psi curves)
% =========================================================
psi_f_abs = abs(psi_f);
psi_b_abs = abs(psi_b);

[psi_f_abs_plot, psi_b_abs_plot, gamma_used] = apply_lorentz_on_doping_( ...
    dop_fwd, psi_f_abs, dop_bwd, psi_b_abs, par.smooth);

% chi plot uses chi_art raw points (dop_use, chi_use)
chi_plot_x = dop_use;
chi_plot_y = chi_use;

% =========================================================
% 7) plot
% =========================================================
fs = par.plot.fontSize;

fig = figure('Color','w','Units','pixels','Position',[80 80 1200 720], ...
    'Name','Hysteresis (shift-averaged chi)');
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

tShort = sprintf('$T=%.6g$, $\\mathrm{artificial\\_polar}=%.6g\\,\\mathrm{meV}$, $q=(%d,%d)$', ...
    T_target, pmev, iq_pick, jq_pick);

% chi(dop) (shift-averaged)
ax1 = nexttile(tl,1);
plot(ax1, chi_plot_x, chi_plot_y, 'o-','LineWidth',2,'MarkerSize',6);
grid(ax1,'on'); box(ax1,'on');
xlabel(ax1,'$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$','Interpreter','latex');
ylabel(ax1,'$\chi_{\mathrm{art}}=\mathrm{Re}\,\chi$','Interpreter','latex');
title(ax1, tShort, 'Interpreter','latex','FontWeight','normal','FontSize',fs);
set(ax1,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out');

% a1(dop) from chi_art
ax2 = nexttile(tl,3);
plot(ax2, dop_fwd, a1_f, 'k-','LineWidth',2);
grid(ax2,'on'); box(ax2,'on');
yline(ax2,0,'--','LineWidth',1.2);
xlabel(ax2,'$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$','Interpreter','latex');
ylabel(ax2,'$a_1(\chi_{\mathrm{art}})$','Interpreter','latex');
set(ax2,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out');

% |psi|
ax3 = nexttile(tl,2);
plot(ax3, dop_fwd, psi_f_abs_plot(:,1), '-','LineWidth',2); hold(ax3,'on');
plot(ax3, dop_bwd, psi_b_abs_plot(:,1), '--','LineWidth',2); hold(ax3,'off');
grid(ax3,'on'); box(ax3,'on');
xlabel(ax3,'$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$','Interpreter','latex');
ylabel(ax3,'$|\psi|$','Interpreter','latex');
title(ax3,'$|\psi|$ hysteresis','Interpreter','latex','FontWeight','normal','FontSize',fs);
legend(ax3,{'forward','backward'},'Interpreter','latex','Location','best');
set(ax3,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out');

% |X|
ax4 = nexttile(tl,4);
plot(ax4, dop_fwd, psi_f_abs_plot(:,2), '-','LineWidth',2); hold(ax4,'on');
plot(ax4, dop_bwd, psi_b_abs_plot(:,2), '--','LineWidth',2); hold(ax4,'off');
grid(ax4,'on'); box(ax4,'on');
xlabel(ax4,'$\mathrm{doping}\ (10^{12}\ \mathrm{cm}^{-2})$','Interpreter','latex');
ylabel(ax4,'$|X|$','Interpreter','latex');
title(ax4,'$|X|$ hysteresis','Interpreter','latex','FontWeight','normal','FontSize',fs);
legend(ax4,{'forward','backward'},'Interpreter','latex','Location','best');
set(ax4,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out');

linkaxes([ax1 ax2 ax3 ax4],'x');

% save
out_dir = fullfile(root, "..", par.plot.save_dir_name);
if ~exist(out_dir,"dir"), mkdir(out_dir); end

out_png = fullfile(out_dir, sprintf("hyst_shiftAvg_T%.4f_artPolar%.6gmeV_iq%d_jq%d_lambda%.3g.png", ...
    T_target, pmev, iq_pick, jq_pick, par.lambda));
exportgraphics(fig, out_png, 'Resolution', par.plot.export_dpi);
fprintf("Saved: %s\n", out_png);

% output
H = struct();
H.root = root;
H.T = T_target;
H.iq_pick = iq_pick;
H.jq_pick = jq_pick;
H.artificial_polar = pmev;
H.dmu_eV = dmu;

H.mu_raw  = mu0;
H.dop_raw = dop0;
H.chi_raw = chi0;

H.mu_used  = mu_use;
H.dop_used = dop_use;
H.chi_used = chi_use;

H.dop_fwd = dop_fwd;
H.dop_bwd = dop_bwd;

H.psi_fwd = psi_f;
H.psi_bwd = psi_b;

H.a1_fwd  = a1_f;
H.a1_bwd  = a1_b;

H.chi_fwd = chi_f;
H.chi_bwd = chi_b;

H.png = out_png;
H.par = par;
end

% =====================================================================
% Collect (mu, doping, chi) at ONE T and ONE q from many chi*.txt
% We accept T from folder name first; header T as fallback.
% mu is parsed from folder name "mu{...}" (preferred) or header.
% doping is parsed from header (required).
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

    % T from path, else header
    T_path = parse_T_from_path_(folder);
    H = parse_header_mu_dop_T_(fpath);

    if isfinite(T_path)
        T_use = T_path;
    else
        T_use = H.T;
    end
    if ~isfinite(T_use), dropH=dropH+1; continue; end
    if ~(T_use>=Tlo && T_use<=Thi), dropT=dropT+1; continue; end

    % mu from path, else header
    mu_path = parse_mu_from_path_(folder);
    mu_use = mu_path;
    if ~isfinite(mu_use), mu_use = H.mu; end

    if ~isfinite(mu_use) || ~isfinite(H.dop)
        dropH = dropH + 1;
        continue;
    end

    % read numeric block
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
tok = regexp(folder, "/T\s*=?\s*" + num + "(?:/|$)", "tokens","once","ignorecase");
if ~isempty(tok)
    T = str2double(tok{1});
else
    tok = regexp(folder, "/temp\s*=?\s*" + num + "(?:/|$)", "tokens","once","ignorecase");
    if ~isempty(tok), T = str2double(tok{1}); else, T = NaN; end
end
end

function mu = parse_mu_from_path_(folder)
folder = replace(string(folder),"\","/");
num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";
tok = regexp(folder, "/mu\s*=?\s*" + num + "(?:/|$)", "tokens","once","ignorecase");
if ~isempty(tok), mu = str2double(tok{1}); else, mu = NaN; end
end

function C = detect_cols_(M)
% supports your two formats
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
% Build hysteresis doping grid from available doping points.
% Supports:
%   hyst.grid_mode = "N"     with hyst.N
%   hyst.grid_mode = "step"  with hyst.step
%
% New:
%   hyst.forward_direction = "ascend"  -> forward: min(dop) -> max(dop)
%   hyst.forward_direction = "descend" -> forward: max(dop) -> min(dop)

dop0 = double(dop0(:));
dop0 = dop0(isfinite(dop0));

if isempty(dop0)
    error("build_hyst_doping_grid_: dop0 is empty.");
end

dmin = min(dop0);
dmax = max(dop0);

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

% -------- build base ascending grid first --------
if mode == "step"
    step = double(hyst.step);
    if ~isfinite(step) || step <= 0
        error("hyst.step must be > 0 when grid_mode='step'");
    end

    dop_base = (dmin:step:dmax).';
    if isempty(dop_base) || abs(dop_base(end) - dmax) > 1e-12
        dop_base(end+1,1) = dmax; %#ok<AGROW>
    end

elseif mode == "N"
    N = double(hyst.N);
    if ~isfinite(N) || N < 3
        error("hyst.N must be >= 3 when grid_mode='N'");
    end

    dop_base = linspace(dmin, dmax, round(N)).';

else
    error("Unsupported hyst.grid_mode. Use 'N' or 'step'.");
end

% -------- apply forward direction --------
if direction == "ascend"
    dop_fwd = dop_base;
else
    dop_fwd = flipud(dop_base);
end

% backward is always reverse of forward
dop_bwd = flipud(dop_fwd);
end

% =========================
% REQUIRED HELPERS (MIN)
% =========================

function [Yf_plot, Yb_plot, gamma_used] = apply_lorentz_on_doping_(xf, Yf, xb, Yb, smooth)
% Plotting-only Lorentz smoothing along x (doping).
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
    return; % OFF
end

if gfac == 0
    gamma = 1.5 * dx; % auto
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
% Continuation: pick the minimum closest to prev
d = sqrt(sum((mins - prev).^2, 2));
[~,ix] = min(d);
pick = mins(ix,:);
end


function mins = find_local_minima_2d_denseSeeds_(F, psi1_lim, psi2_lim, scan, opt)
% Dense scan -> local-min mask -> pick seeds -> fminsearch -> cluster -> Hessian PD

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