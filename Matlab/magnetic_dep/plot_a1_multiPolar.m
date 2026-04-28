function H = plot_a1_multiPolar_vs_x()
% Plot a1 curves for multiple artificial polar shifts.
%
% Core rule (same as user's hysteresis code):
%   dmu = 0.5 * artificial_polar(meV) * 1e-3  [eV]
%   chi_art(mu) = 0.5 * [ chi(mu+dmu) + chi(mu-dmu) ]
%
% Then compute:
%   C = coef.eval(T_target, doping, chi_art)
%   a1 = C.a1
%
% Features:
%   - multiple polar values (meV)
%   - x-axis selectable: doping or mu
%   - hot colormap, but flipped so 0 meV = white end
%   - actual 0 meV curve is drawn as gray dashed
%   - recursive read of chi*.txt under chosen root
%   - no auto save

% ---------------- params / coeff ----------------
par  = make_realchi_params(true);
coef = make_realchi_coeff(par);

% ---------------- user settings -----------------
default_root = "E:/rg_master/data/";
if ~isfolder(default_root)
    default_root = string(pwd);
end

polar_list_meV = [0];   % <--- edit here
FS = 14;
LW0 = 2.2;
LW  = 2.2;

% ---------------- choose x-axis -----------------
x_mode = pick_xmode_();   % "doping" or "mu"
fprintf("[x_mode] %s\n", x_mode);

% ---------------- choose folder -----------------
root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
if isequal(root, 0), error('User cancelled.'); end
root = string(root);
fprintf("Root: %s\n", root);

% ---------------- q point -----------------------
iq_pick = round(par.iq_pick);
jq_pick = round(par.jq_pick);
fprintf("Pick (iq,jq)=(%d,%d)\n", iq_pick, jq_pick);

% ---------------- target T ----------------------
T_target = par.T_target;
T_tol    = par.T_tol;
fprintf("T_target=%.6g K, T_tol=%.3g K\n", T_target, T_tol);

% ---------------- collect baseline --------------
B = collect_mu_dop_chi_curve_(root, T_target, T_tol, iq_pick, jq_pick);

mu0  = B.mu(:);
dop0 = B.dop(:);
chi0 = B.chi(:);

if isempty(mu0)
    fprintf("[baseline] points=0 | mu=[, ] | dop=[, ]\n");
else
    fprintf("[baseline] points=%d | mu=[%.6g, %.6g] | dop=[%.6g, %.6g]\n", ...
        numel(mu0), min(mu0), max(mu0), min(dop0), max(dop0));
end

if numel(mu0) < 5
    error("Too few points at this T. Check T_tol / header parse / (iq,jq).");
end

% sort by mu & collapse duplicates
[mu0, chi0, dop0] = collapse_same_mu_(mu0, chi0, dop0);

% build interpolants on mu
Fchi0 = griddedInterpolant(mu0, chi0, "linear", "none");
Fdop0 = griddedInterpolant(mu0, dop0, "linear", "none");

% ---------------- sort polar list ---------------
polar_list_meV = unique(polar_list_meV(:).', 'sorted');
fprintf("[polar list meV] ");
fprintf("%.6g ", polar_list_meV);
fprintf("\n");

% ---------------- colormap ----------------------
% hot default: black -> red -> yellow -> white
% flip it so 0 meV corresponds to white end
pmin = min(polar_list_meV);
pmax = max(polar_list_meV);

n_cmap = 256;
cmap_full = flipud(hot(n_cmap));

if pmin == pmax
    get_polar_color = @(p) cmap_full(1,:);
else
    get_polar_color = @(p) interp1( ...
        linspace(pmin, pmax, n_cmap), ...
        cmap_full, ...
        p, ...
        'linear', 'extrap');
end

% ---------------- figure ------------------------
fig = figure('Color','w','Units','pixels','Position',[120 120 860 580], ...
    'Name', sprintf('a1 vs %s | multi-polar', x_mode));

ax = axes(fig);
hold(ax,'on');
box(ax,'on');
grid(ax,'on');
set(ax,'FontSize',FS,'LineWidth',1,'TickDir','out');

% store outputs
H = struct();
H.root = root;
H.T_target = T_target;
H.T_tol = T_tol;
H.iq_pick = iq_pick;
H.jq_pick = jq_pick;
H.x_mode = x_mode;

H.mu_raw = mu0;
H.dop_raw = dop0;
H.chi_raw = chi0;

curve_data = struct([]);

% ---------------- loop over polar ---------------
for ip = 1:numel(polar_list_meV)
    pmev = polar_list_meV(ip);
    dmu  = 0.5 * pmev * 1e-3;   % eV

    fprintf("[polar] %.6g meV => dmu = %.12g eV\n", pmev, dmu);

    % ---- build chi_art(mu) on original mu grid
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
        fprintf("[skip] polar=%.6g meV | too few finite points after shift-average\n", pmev);
        continue;
    end

    % ---- compute a1 pointwise
    a1_use = nan(size(mu_use));
    for i = 1:numel(mu_use)
        C = coef.eval(T_target, dop_use(i), chi_use(i));
        a1_use(i) = C.a1;
    end

    % ---- choose x-axis
    if x_mode == "doping"
        x_use = dop_use;
        [x_plot, a1_plot, mu_plot, chi_plot] = collapse_same_x_(x_use, a1_use, mu_use, chi_use);
        [x_plot, ord] = sort(x_plot);
        a1_plot = a1_plot(ord);
        mu_plot = mu_plot(ord);
        chi_plot = chi_plot(ord);
        x_label = 'doping (10^{12} cm^{-2})';
    else
        x_use = mu_use;
        [x_plot, a1_plot, dop_plot, chi_plot] = collapse_same_x_(x_use, a1_use, dop_use, chi_use);
        [x_plot, ord] = sort(x_plot);
        a1_plot = a1_plot(ord);
        dop_plot = dop_plot(ord);
        chi_plot = chi_plot(ord);
        x_label = '\mu (eV)';
    end

    % ---- plot
    if pmev == 0
        plot(ax, x_plot, a1_plot, '--', ...
            'Color', [0.45 0.45 0.45], ...
            'LineWidth', LW0, ...
            'DisplayName', 'a_1 (polar = 0 meV)');
    else
        this_color = get_polar_color(pmev);
        plot(ax, x_plot, a1_plot, '-', ...
            'Color', this_color, ...
            'LineWidth', LW, ...
            'DisplayName', sprintf('polar = %.4g meV', pmev));
    end

    % ---- save curve data
    curve_data(end+1).polar_meV = pmev; %#ok<AGROW>
    curve_data(end).dmu_eV = dmu;
    curve_data(end).x = x_plot;
    curve_data(end).a1 = a1_plot;
    if x_mode == "doping"
        curve_data(end).mu = mu_plot;
    else
        curve_data(end).dop = dop_plot;
    end
    curve_data(end).chi = chi_plot;
end

% ---------------- labels ------------------------
xlabel(ax, x_label, 'Interpreter','tex', 'FontSize',FS);
ylabel(ax, '$a_1$', 'Interpreter','latex', 'FontSize',FS);

tShort = sprintf('$T=%.6g\\,\\mathrm{K}$, $q=(%d,%d)$', T_target, iq_pick, jq_pick);
title(ax, ['$a_1$ for multiple polar shifts, ' tShort], ...
    'Interpreter','latex', 'FontWeight','normal', 'FontSize', FS);

yline(ax, 0, '--', 'LineWidth', 1.2);

% ---------------- colorbar ----------------------
colormap(ax, cmap_full);
cb = colorbar(ax);
cb.Label.String = 'polar (meV)';
cb.FontSize = FS;

if pmin == pmax
    caxis(ax, [pmin-1, pmax+1]);
else
    caxis(ax, [pmin, pmax]);
end
cb.Ticks = polar_list_meV;

legend(ax, 'Location', 'best');
set(ax,'TickLabelInterpreter','latex');

% ---------------- output ------------------------
H.fig = fig;
H.ax = ax;
H.curves = curve_data;
H.polar_list_meV = polar_list_meV;
H.par = par;

end

% ============================================================
% helper: choose x-axis
% ============================================================
function x_mode = pick_xmode_()
c = questdlg('Choose x-axis:', 'x-axis', 'doping', 'mu', 'doping');
if isempty(c)
    x_mode = "doping";
else
    x_mode = string(c);
end
end

% ============================================================
% collect (mu, doping, chi) at ONE T and ONE q from chi*.txt
% ============================================================
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

function [x2, y2, z2, w2] = collapse_same_x_(x, y, z, w)
x = double(x(:));
y = double(y(:));
z = double(z(:));
w = double(w(:));

ok = isfinite(x) & isfinite(y) & isfinite(z) & isfinite(w);
x=x(ok); y=y(ok); z=z(ok); w=w(ok);

[xu,~,ic] = unique(x);
ym = accumarray(ic, y, [], @mean);
zm = accumarray(ic, z, [], @mean);
wm = accumarray(ic, w, [], @mean);

[x2,ord] = sort(xu);
y2 = ym(ord);
z2 = zm(ord);
w2 = wm(ord);
end