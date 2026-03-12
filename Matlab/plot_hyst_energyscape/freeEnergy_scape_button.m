function freeEnergy_scape_button()
% Minimal standalone UI: pick folder -> parse T & (mu/doping) from path -> read chi -> plot F(psi1,psi2)

par  = make_realchi_params();
coef = make_realchi_coeff(par);

% choose folder
default_root = string(pwd);
if isfield(par,'io') && isfield(par.io,'default_root') && isfolder(par.io.default_root)
    default_root = string(par.io.default_root);
end
root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
if isequal(root,0), return; end
root = string(root);

% load grid (embedded)
G = load_chi_grid_folderUI_(root, par.iq_pick, par.jq_pick);

T_list  = G.T_list;
U_list  = G.U_list;
chi_map = G.chi_map;
dop_map = G.doping_map;
u_tag   = G.u_tag;

% style
FS  = 14; dpi = 300;
if isfield(par,'plot')
    if isfield(par.plot,'fontSize'), FS = par.plot.fontSize; end
    if isfield(par.plot,'export_dpi'), dpi = par.plot.export_dpi; end %#ok<NASGU>
end
Ngrid2D = 161;

% UI
fig = uifigure('Name', sprintf('Free Energy (UI=%s)', u_tag), 'Position',[80 80 980 600]);

ax2d = uiaxes(fig,'Position',[50 200 520 400],'FontSize',FS);
ax2d.TickLabelInterpreter = 'latex';
xlabel(ax2d,'$\psi_{\rm CDW}$','Interpreter','latex');
ylabel(ax2d,'$\psi_{\rm lattice}$','Interpreter','latex');
xlim(ax2d, par.psi1_lim); ylim(ax2d, par.psi2_lim);

axPsi = uiaxes(fig,'Position',[600 410 340 190],'FontSize',FS);
axPsi.TickLabelInterpreter = 'latex';
xlabel(axPsi,'$\psi_1$','Interpreter','latex');
ylabel(axPsi,'$F_\psi(\psi_1)$','Interpreter','latex');
title(axPsi,'$F_\psi=\frac12 a_1\psi^2+\frac{1}{4!}b_1\psi^4$','Interpreter','latex','FontWeight','normal');

axX = uiaxes(fig,'Position',[600 200 340 190],'FontSize',FS);
axX.TickLabelInterpreter = 'latex';
xlabel(axX,'$\psi_2$','Interpreter','latex');
ylabel(axX,'$F_X(\psi_2)$','Interpreter','latex');
title(axX,'$F_X=\frac12 a_2X^2+\frac{1}{4!}b_2X^4+\frac{1}{6!}c_2X^6$','Interpreter','latex','FontWeight','normal');

boxC = uitextarea(fig,'Position',[600 40 340 140],'Editable','off','FontSize',FS-2);

% state
[it, iu] = find_first_valid_(chi_map);
lambda = par.lambda0;
psi1_opt = 0; psi2_opt = 0;

% controls
x0=70; y1=95; dy=48; btnW=34; btnH=28; gapB=6; gap=12; valW=320; valH=24;
mk4 = @(y) deal( ...
    uibutton(fig,'push','Position',[x0+50+0*(btnW+gapB) y btnW btnH],'Text','--'), ...
    uibutton(fig,'push','Position',[x0+50+1*(btnW+gapB) y btnW btnH],'Text','-'),  ...
    uibutton(fig,'push','Position',[x0+50+2*(btnW+gapB) y btnW btnH],'Text','+'),  ...
    uibutton(fig,'push','Position',[x0+50+3*(btnW+gapB) y btnW btnH],'Text','++') ...
);

uilabel(fig,'Position',[x0 y1+2*dy 60 22],'Text','T','FontSize',FS-2);
[Tmm,Tm,Tp,Tpp] = mk4(y1+2*dy-4);
labT = uilabel(fig,'Position',[x0+50+4*(btnW+gapB)+gap y1+2*dy-2 valW valH],'Text','','FontSize',FS-2,'HorizontalAlignment','left');

uilabel(fig,'Position',[x0 y1+dy 80 22],'Text',u_tag,'FontSize',FS-2);
[Umm,Um,Up,Upp] = mk4(y1+dy-4);
labU = uilabel(fig,'Position',[x0+50+4*(btnW+gapB)+gap y1+dy-2 valW valH],'Text','','FontSize',FS-2,'HorizontalAlignment','left');

uilabel(fig,'Position',[x0 y1 60 22],'Text','lambda','FontSize',FS-2);
[Lmm,Lm,Lp,Lpp] = mk4(y1-4);
labL = uilabel(fig,'Position',[x0+50+4*(btnW+gapB)+gap y1-2 valW valH],'Text','','FontSize',FS-2,'HorizontalAlignment','left');

iq_in = uieditfield(fig,'numeric','Position',[200 20 70 30],'Value',par.iq_pick,'FontSize',FS-2);
jq_in = uieditfield(fig,'numeric','Position',[280 20 70 30],'Value',par.jq_pick,'FontSize',FS-2);
uilabel(fig,'Position',[200 50 70 22],'Text','iq','FontSize',FS-2);
uilabel(fig,'Position',[280 50 70 22],'Text','jq','FontSize',FS-2);
btnApply = uibutton(fig,'push','Position',[360 20 110 30],'Text','Apply iq/jq','FontSize',FS-2);

bigT = max(3, round(numel(T_list)/20));
bigU = 10;
dLamS = 0.1; dLamB = 0.5;

Tmm.ButtonPushedFcn = @(~,~) stepT(-bigT);
Tm.ButtonPushedFcn  = @(~,~) stepT(-1);
Tp.ButtonPushedFcn  = @(~,~) stepT(+1);
Tpp.ButtonPushedFcn = @(~,~) stepT(+bigT);

Umm.ButtonPushedFcn = @(~,~) stepU(-bigU);
Um.ButtonPushedFcn  = @(~,~) stepU(-1);
Up.ButtonPushedFcn  = @(~,~) stepU(+1);
Upp.ButtonPushedFcn = @(~,~) stepU(+bigU);

Lmm.ButtonPushedFcn = @(~,~) stepLam(-dLamB);
Lm.ButtonPushedFcn  = @(~,~) stepLam(-dLamS);
Lp.ButtonPushedFcn  = @(~,~) stepLam(+dLamS);
Lpp.ButtonPushedFcn = @(~,~) stepLam(+dLamB);

btnApply.ButtonPushedFcn = @(~,~) applyIQJQ();

refresh(); drawAll();

% ---------------- nested ----------------
    function stepT(d)
        it = clamp_(it + d, numel(T_list));
        iu = nearest_valid_in_row_(chi_map, it, iu);
        refresh(); drawAll();
    end

    function stepU(d)
        iu = clamp_(iu + d, numel(U_list));
        iu = nearest_valid_in_row_(chi_map, it, iu);
        refresh(); drawAll();
    end

    function stepLam(d)
        lambda = lambda + d;
        lambda = round(lambda/dLamS)*dLamS;
        refresh(); drawAll();
    end

    function refresh()
        labT.Text = sprintf('%.6g K', T_list(it));
        if strcmpi(u_tag,'doping'), labU.Text = sprintf('%.4f', U_list(iu));
        else, labU.Text = sprintf('%.6g', U_list(iu));
        end
        labL.Text = sprintf('%.2f', lambda);
    end

    function drawAll()
        if ~isfinite(chi_map(it,iu))
            iu = nearest_valid_in_row_(chi_map, it, iu);
        end

        T = T_list(it);
        u = U_list(iu);
        chi_used = chi_map(it,iu);
        dop_eval = dop_map(it,iu);

        C = coef.eval(T, dop_eval, chi_used);

        try
            [psi1_opt, psi2_opt] = minimize_free_energy( ...
                C.a1, C.b1, C.a2, C.b2, C.c2, lambda, ...
                psi1_opt + 1e-3*randn(), psi2_opt + 1e-3*randn(), false);
        catch
            psi1_opt = 0; psi2_opt = 0;
        end

        F2D = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, lambda);

        cla(ax2d); cla(axPsi); cla(axX);

        if strcmpi(u_tag,'doping'), ustr = sprintf('%.4f', u);
        else, ustr = sprintf('%.6g', u);
        end
        title(ax2d, sprintf('T=%.6g, %s=%s, lambda=%.2f', T, u_tag, ustr, lambda), ...
            'Interpreter','none','FontWeight','normal');

        psi1_vals = linspace(par.psi1_lim(1), par.psi1_lim(2), Ngrid2D);
        psi2_vals = linspace(par.psi2_lim(1), par.psi2_lim(2), Ngrid2D);
        [P1,P2] = meshgrid(psi1_vals, psi2_vals);
        Fv = arrayfun(F2D, P1, P2);

        pcolor(ax2d, P1, P2, Fv); shading(ax2d,'interp'); colormap(ax2d,'turbo');
        colorbar(ax2d,'Location','eastoutside','FontSize',FS-2);
        hold(ax2d,'on');
        contour(ax2d, P1, P2, Fv, 22, 'LineColor','k','LineWidth',0.5);
        plot(ax2d, psi1_opt, psi2_opt, 'ro','MarkerFaceColor','r','MarkerSize',7);
        hold(ax2d,'off');
        xlim(ax2d, par.psi1_lim); ylim(ax2d, par.psi2_lim);

        Fpsi = @(psi) 0.5*C.a1*psi.^2 + (1/factorial(4))*C.b1*psi.^4;
        FX   = @(X)   0.5*C.a2*X.^2   + (1/factorial(4))*C.b2*X.^4 + (1/factorial(6))*C.c2*X.^6;

        x1 = linspace(par.psi1_lim(1), par.psi1_lim(2), 801);
        x2 = linspace(par.psi2_lim(1), par.psi2_lim(2), 801);

        plot(axPsi, x1, Fpsi(x1), 'LineWidth',2); grid(axPsi,'on'); box(axPsi,'on');
        hold(axPsi,'on'); plot(axPsi, psi1_opt, Fpsi(psi1_opt), 'ro','MarkerFaceColor','r','MarkerSize',6); hold(axPsi,'off');

        plot(axX, x2, FX(x2), 'LineWidth',2); grid(axX,'on'); box(axX,'on');
        hold(axX,'on'); plot(axX, psi2_opt, FX(psi2_opt), 'ro','MarkerFaceColor','r','MarkerSize',6); hold(axX,'off');

        boxC.Value = {
            sprintf('T = %.12g', T)
            sprintf('%s(folder) = %.12g', u_tag, u)
            sprintf('doping(header) = %.12g', dop_eval)
            sprintf('chi = %.12g', chi_used)
            sprintf('a1=%.6g  b1=%.6g', C.a1, C.b1)
            sprintf('a2=%.6g  b2=%.6g  c2=%.6g', C.a2, C.b2, C.c2)
            sprintf('psi1*=%.6g  psi2*=%.6g', psi1_opt, psi2_opt)
        };
    end

    function applyIQJQ()
        par.iq_pick = round(iq_in.Value);
        par.jq_pick = round(jq_in.Value);

        G = load_chi_grid_folderUI_(root, par.iq_pick, par.jq_pick);
        T_list  = G.T_list;
        U_list  = G.U_list;
        chi_map = G.chi_map;
        dop_map = G.doping_map;
        u_tag   = G.u_tag;

        [it, iu] = find_first_valid_(chi_map);
        psi1_opt = 0; psi2_opt = 0;
        refresh(); drawAll();
    end
end

% =====================================================================
% Embedded loader: parse folder Txxx & muXXX/dopingXXX, read header T+doping,
% read numeric, pick Re chi(iq,jq), average duplicates.
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

u_tag = "";

for k = 1:n
    fpath = string(fullfile(L(k).folder, L(k).name));

    [~,~,ukind,uval] = parse_TU_from_path_(fpath);
    if ukind=="" || ~isfinite(uval), continue; end
    if u_tag=="" && ukind~="", u_tag = ukind; end

    H = parse_header_Tdop_(fpath);
    if ~H.ok, continue; end

    M = read_numeric_skiphash_(fpath);
    if isempty(M) || size(M,2) < 5, continue; end

    Cc = detect_cols_(M);
    id = find(M(:,Cc.iq)==iq_pick & M(:,Cc.jq)==jq_pick, 1);
    if isempty(id), continue; end

    Tvals(k) = round(H.T_K, 12);
    Uvals(k) = quantize_U_(ukind, uval);
    chiV(k)  = M(id, Cc.Re);
    dopH(k)  = round(H.doping, 12);
end

if u_tag=="", u_tag="doping"; end

mask = isfinite(Tvals) & isfinite(Uvals) & isfinite(chiV) & isfinite(dopH);
Tvals=Tvals(mask); Uvals=Uvals(mask); chiV=chiV(mask); dopH=dopH(mask);

T_list = sort(unique(Tvals));
U_list = sort(unique(Uvals));
NT = numel(T_list); NU = numel(U_list);

chi_sum = zeros(NT,NU); chi_cnt = zeros(NT,NU);
dop_sum = zeros(NT,NU); dop_cnt = zeros(NT,NU);

for i = 1:numel(Tvals)
    iT = find(T_list==Tvals(i),1);
    iU = find(U_list==Uvals(i),1);
    if isempty(iT)||isempty(iU), continue; end
    chi_sum(iT,iU)=chi_sum(iT,iU)+chiV(i); chi_cnt(iT,iU)=chi_cnt(iT,iU)+1;
    dop_sum(iT,iU)=dop_sum(iT,iU)+dopH(i); dop_cnt(iT,iU)=dop_cnt(iT,iU)+1;
end

chi_map = nan(NT,NU); m=chi_cnt>0; chi_map(m)=chi_sum(m)./chi_cnt(m);
doping_map = nan(NT,NU); md=dop_cnt>0; doping_map(md)=dop_sum(md)./dop_cnt(md);

G = struct('T_list',T_list,'U_list',U_list,'chi_map',chi_map,'doping_map',doping_map,'u_tag',u_tag);
end

function [Tfolder, Ufolder, ukind, uval] = parse_TU_from_path_(fpath)
Tfolder=""; Ufolder=""; ukind=""; uval=NaN;
p = replace(string(fpath), "\", "/");
tokT = regexp(p, '/(T[^/]+)/', 'tokens', 'once'); if ~isempty(tokT), Tfolder=string(tokT{1}); end %#ok<NASGU>
extract_num = @(s) str2double(regexp(s,'([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)','match','once'));

tokMu = regexp(p, '/(mu[^/]+)/', 'tokens', 'once');
if ~isempty(tokMu)
    Ufolder = string(tokMu{1}); ukind="mu"; uval=extract_num(Ufolder);
    if isfinite(uval), return; end
end

tokDp = regexp(p, '/(doping[^/]+)/', 'tokens', 'once');
if ~isempty(tokDp)
    Ufolder = string(tokDp{1}); ukind="doping"; uval=extract_num(Ufolder);
end
end

function u = quantize_U_(ukind, uval)
if ukind=="doping", u = round(uval,4);
else, u = round(uval,12);
end
end

function H = parse_header_Tdop_(fpath)
H = struct('ok',false,'T_K',NaN,'doping',NaN);
fid = fopen(fpath,'r'); if fid<0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';
for t = 1:3000
    ln = fgetl(fid); if ~ischar(ln), break; end
    s = strtrim(ln); if isempty(s), continue; end
    if s(1)~='#', break; end
    if ~isfinite(H.T_K)
        tok = regexp(s, ['^\#\s*T\s*=\s*' num], 'tokens','once');
        if isempty(tok), tok = regexp(s, ['^\#\s*T_K\s*=\s*' num], 'tokens','once'); end
        if ~isempty(tok), H.T_K = str2double(tok{1}); end
    end
    if ~isfinite(H.doping)
        tok = regexp(s, ['^\#\s*doping\s*=\s*' num], 'tokens','once');
        if ~isempty(tok), H.doping = str2double(tok{1}); end
    end
    if isfinite(H.T_K) && isfinite(H.doping), break; end
end
H.ok = isfinite(H.T_K) && isfinite(H.doping);
end

function M = read_numeric_skiphash_(fpath)
fid = fopen(fpath,'r'); if fid<0, M=[]; return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
rows = {};
while true
    ln = fgetl(fid); if ~ischar(ln), break; end
    if isempty(ln), continue; end
    if ~isempty(regexp(ln,'^\s*#','once')), continue; end
    v = sscanf(ln,'%f').';
    if isempty(v), continue; end
    rows{end+1,1} = v; %#ok<AGROW>
end
if isempty(rows), M=[]; return; end
ncol = max(cellfun(@numel, rows));
M = nan(numel(rows), ncol);
for i=1:numel(rows)
    v=rows{i}; M(i,1:numel(v))=v;
end
lastFinite = find(any(isfinite(M),1),1,'last');
if ~isempty(lastFinite), M=M(:,1:lastFinite); end
end

function C = detect_cols_(M)
ncol = size(M,2);
C = struct('iq',1,'jq',2,'Re',5);
if ncol >= 8
    C.iq = 2; C.jq = 3; C.Re = 6;
end
end

% ====================== tiny helpers ======================
function i = clamp_(i, N)
i = max(1, min(N, i));
end

function [it, iu] = find_first_valid_(M)
[it, iu] = find(isfinite(M), 1, 'first');
if isempty(it), it=1; iu=1; end
end

function iu2 = nearest_valid_in_row_(M, it, iu0)
row = M(it,:);
good = find(isfinite(row));
if isempty(good), iu2 = iu0; return; end
[~,k] = min(abs(good - iu0));
iu2 = good(k);
end