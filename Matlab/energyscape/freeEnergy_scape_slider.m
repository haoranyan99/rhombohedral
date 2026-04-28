function freeEnergy_scape_slider()
% Minimal standalone SLIDER UI version (single file).
% Depends on your existing:
%   make_realchi_params, make_realchi_coeff, free_energy, minimize_free_energy

par  = make_realchi_params(true);
coef = make_realchi_coeff(par);

% folder
% default_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_200\D0.084";
default_root = "/Users/haoranyan/rg_master/data/";
if isfield(par,'io') && isfield(par.io,'default_root') && isfolder(par.io.default_root)
    default_root = string(par.io.default_root);
end
root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
if isequal(root,0), return; end
root = string(root);

% load (embedded)
G = load_chi_grid_folderUI_(root, par.iq_pick, par.jq_pick);
T_list  = G.T_list;
U_list  = G.U_list;
chi_map = G.chi_map;
dop_map = G.doping_map;
muf_map = G.mu_folder_map;
u_tag   = G.u_tag;

[it0, iu0] = find_first_valid_(chi_map);
psi1_opt = 0; psi2_opt = 0;

% UI
fig = uifigure('Name', sprintf('Free Energy (slider UI=%s)', u_tag), 'Position',[80 80 980 600]);

ax2d = uiaxes(fig,'Position',[50 200 520 400],'FontSize',12);
ax2d.TickLabelInterpreter = 'latex';
xlabel(ax2d,'$\psi_{\rm CDW}$','Interpreter','latex');
ylabel(ax2d,'$\psi_{\rm lattice}$','Interpreter','latex');
xlim(ax2d, par.psi1_lim); ylim(ax2d, par.psi2_lim);

axPsi = uiaxes(fig,'Position',[600 410 340 190],'FontSize',12);
axPsi.TickLabelInterpreter = 'latex';
xlabel(axPsi,'$\psi_1$','Interpreter','latex');
ylabel(axPsi,'$F_\psi(\psi_1)$','Interpreter','latex');
title(axPsi,'$F_\psi=\frac12 a_1\psi^2+\frac{1}{4!}b_1\psi^4$','Interpreter','latex','FontWeight','normal');

axX = uiaxes(fig,'Position',[600 200 340 190],'FontSize',12);
axX.TickLabelInterpreter = 'latex';
xlabel(axX,'$\psi_2$','Interpreter','latex');
ylabel(axX,'$F_X(\psi_2)$','Interpreter','latex');
title(axX,'$F_X=\frac12 a_2X^2+\frac{1}{4!}b_2X^4+\frac{1}{6!}c_2X^6$','Interpreter','latex','FontWeight','normal');

coeff_box = uitextarea(fig,'Position',[600 40 340 140],'Editable','off','FontSize',12);

% sliders
y1=90; dy=50;

T_slider = uislider(fig,'Position',[100 y1+2*dy 330 3], ...
    'Limits',[1, max(1,numel(T_list))], 'Value', it0);
T_slider.MajorTicks=[]; T_slider.MinorTicks=[];

U_slider = uislider(fig,'Position',[100 y1+dy 330 3], ...
    'Limits',[1, max(1,numel(U_list))], 'Value', iu0);
U_slider.MajorTicks=[]; U_slider.MinorTicks=[];

% --- lambda slider range: [0, par.lambda] ---
lam_max = par.lambda;
if isempty(lam_max) || ~isfinite(lam_max) || lam_max <= 0
    lam_max = 1; % fallback
end

lambda0 = par.lambda0;
if isempty(lambda0) || ~isfinite(lambda0)
    lambda0 = 0;
end
lambda0 = min(lam_max, max(0, lambda0));

lambda_slider = uislider(fig,'Position',[100 y1 330 3], ...
    'Limits',[0, lam_max], 'Value', lambda0);
lambda_slider.MajorTicks = linspace(0, lam_max, 6);
lambda_slider.MinorTicks = [];
uilabel(fig,'Position',[70 y1+2*dy-10 50 22],'Text','$T$','FontSize',12,'Interpreter','latex');
uilabel(fig,'Position',[70 y1+dy-10  70 22],'Text',sprintf('$%s$',u_tag),'FontSize',12,'Interpreter','latex');
uilabel(fig,'Position',[70 y1-10     70 22],'Text','$\lambda$','FontSize',12,'Interpreter','latex');

T_text = uilabel(fig,'Position',[450 y1+2*dy-14 200 22],'Text','', 'FontSize',12);
U_text = uilabel(fig,'Position',[450 y1+dy-14  240 22],'Text','', 'FontSize',12);
L_text = uilabel(fig,'Position',[450 y1-14     200 22],'Text','', 'FontSize',12);

% iq/jq apply
iq_input = uieditfield(fig,'numeric','Position',[180 20 60 30],'Value',par.iq_pick,'FontSize',12);
jq_input = uieditfield(fig,'numeric','Position',[250 20 60 30],'Value',par.jq_pick,'FontSize',12);
uilabel(fig,'Position',[180 50 60 22],'Text','iq','FontSize',12);
uilabel(fig,'Position',[250 50 60 22],'Text','jq','FontSize',12);
apply_button = uibutton(fig,'push','Position',[320 20 110 30],'Text','Apply iq/jq','FontSize',12);

T_slider.ValueChangedFcn      = @(~,~) update_();
U_slider.ValueChangedFcn      = @(~,~) update_();
lambda_slider.ValueChangedFcn = @(~,~) update_();
apply_button.ButtonPushedFcn  = @(~,~) apply_iqjq_();

refresh_labels_();
update_plot_only_(); % no minimization jump at start

% ========================= nested =========================
    function refresh_labels_()
        it = clamp_(round(T_slider.Value), numel(T_list));
        iu = clamp_(round(U_slider.Value), numel(U_list));
        T_slider.Value = it; U_slider.Value = iu;

        T_text.Text = sprintf('T = %.6g K', T_list(it));
        uval = U_list(iu);
        if strcmpi(u_tag,'doping')
            U_text.Text = sprintf('doping(folder) = %.4f', uval);
        else
            U_text.Text = sprintf('\\mu(folder) = %.6g', uval);
        end
        L_text.Text = sprintf('\\lambda = %.3f', lambda_slider.Value);
    end

    function update_()
        refresh_labels_();
        it = clamp_(round(T_slider.Value), numel(T_list));
        iu = clamp_(round(U_slider.Value), numel(U_list));

        if ~isfinite(chi_map(it,iu))
            iu = nearest_valid_in_row_(chi_map, it, iu);
            U_slider.Value = iu;
            refresh_labels_();
            if ~isfinite(chi_map(it,iu))
                render_empty_(sprintf('No valid point at T=%.6g', T_list(it)));
                return;
            end
        end

        T = T_list(it);
        u = U_list(iu);
        chi_used = chi_map(it,iu);
        dop_eval = dop_map(it,iu);
        mu_folder_value = muf_map(it,iu);
        lambda = lambda_slider.Value;

        C = coef.eval(T, dop_eval, chi_used);

        try
            [psi1_opt, psi2_opt] = minimize_free_energy( ...
                C.a1, C.b1, C.a2, C.b2, C.c2, lambda, ...
                psi1_opt + 1e-3*randn(), psi2_opt + 1e-3*randn(), false);
        catch
            psi1_opt = 0; psi2_opt = 0;
        end

        draw_cached_(T,u,chi_used,dop_eval,mu_folder_value,lambda,C);
    end

    function update_plot_only_()
        refresh_labels_();
        it = clamp_(round(T_slider.Value), numel(T_list));
        iu = clamp_(round(U_slider.Value), numel(U_list));

        if ~isfinite(chi_map(it,iu))
            iu = nearest_valid_in_row_(chi_map, it, iu);
            U_slider.Value = iu;
            refresh_labels_();
            if ~isfinite(chi_map(it,iu))
                render_empty_('No valid initial point.');
                return;
            end
        end

        T = T_list(it);
        u = U_list(iu);
        chi_used = chi_map(it,iu);
        dop_eval = dop_map(it,iu);
        mu_folder_value = muf_map(it,iu);
        lambda = lambda_slider.Value;

        C = coef.eval(T, dop_eval, chi_used);
        draw_cached_(T,u,chi_used,dop_eval,mu_folder_value,lambda,C);
    end

    function draw_cached_(T,u,chi_used,dop_eval,mu_folder_value,lambda,C)
        cla(ax2d); cla(axPsi); cla(axX);

        F2D = free_energy(C.a1, C.b1, C.a2, C.b2, C.c2, lambda);
        psi1_vals = linspace(par.psi1_lim(1), par.psi1_lim(2), 161);
        psi2_vals = linspace(par.psi2_lim(1), par.psi2_lim(2), 161);
        [P1,P2] = meshgrid(psi1_vals, psi2_vals);
        Fv = arrayfun(F2D, P1, P2);

        pcolor(ax2d, P1, P2, Fv); shading(ax2d,'interp'); colormap(ax2d,'turbo');
        colorbar(ax2d,'Location','eastoutside','FontSize',12);
        hold(ax2d,'on');
        contour(ax2d, P1, P2, Fv, 22, 'LineColor','k','LineWidth',0.5);
        plot(ax2d, psi1_opt, psi2_opt, 'ro','MarkerFaceColor','r','MarkerSize',7);
        hold(ax2d,'off');
        xlim(ax2d, par.psi1_lim); ylim(ax2d, par.psi2_lim);

        if strcmpi(u_tag,'doping'), ustr = sprintf('doping(folder)=%.4f', u);
        else, ustr = sprintf('\\mu(folder)=%.6g', u);
        end
        title(ax2d, sprintf('T=%.6g K, %s, \\lambda=%.3f', T, ustr, lambda), ...
            'FontWeight','normal');

        Fpsi = @(psi) 0.5*C.a1*psi.^2 + (1/factorial(4))*C.b1*psi.^4;
        FX   = @(X)   0.5*C.a2*X.^2   + (1/factorial(4))*C.b2*X.^4 + (1/factorial(6))*C.c2*X.^6;

        x1 = linspace(par.psi1_lim(1), par.psi1_lim(2), 801);
        x2 = linspace(par.psi2_lim(1), par.psi2_lim(2), 801);

        plot(axPsi, x1, Fpsi(x1), 'LineWidth',2); grid(axPsi,'on'); box(axPsi,'on');
        hold(axPsi,'on'); plot(axPsi, psi1_opt, Fpsi(psi1_opt), 'ro','MarkerFaceColor','r','MarkerSize',6); hold(axPsi,'off');

        plot(axX, x2, FX(x2), 'LineWidth',2); grid(axX,'on'); box(axX,'on');
        hold(axX,'on'); plot(axX, psi2_opt, FX(psi2_opt), 'ro','MarkerFaceColor','r','MarkerSize',6); hold(axX,'off');

        lines = {
            sprintf('iq = %d,   jq = %d', par.iq_pick, par.jq_pick)
            sprintf('T = %.12g   lambda = %.12g', T, lambda)
            sprintf('UI: %s(folder) = %.12g', u_tag, u)
            sprintf('doping(header) = %.12g', dop_eval)
        };
        if strcmpi(u_tag,'mu')
            lines{end+1} = sprintf('mu(folder) = %.12g', mu_folder_value);
        end
        lines = [lines; {
            sprintf('chi = %.12g', chi_used)
            sprintf('a1 = %.12g', C.a1)
            sprintf('b1 = %.12g', C.b1)
            sprintf('a2 = %.12g', C.a2)
            sprintf('b2 = %.12g', C.b2)
            sprintf('c2 = %.12g', C.c2)
            sprintf('psi1* = %.12g    psi2* = %.12g', psi1_opt, psi2_opt)
        }];
        coeff_box.Value = lines;
    end

    function apply_iqjq_()
        par.iq_pick = round(iq_input.Value);
        par.jq_pick = round(jq_input.Value);

        G = load_chi_grid_folderUI_(root, par.iq_pick, par.jq_pick);
        T_list  = G.T_list;
        U_list  = G.U_list;
        chi_map = G.chi_map;
        dop_map = G.doping_map;
        muf_map = G.mu_folder_map;
        u_tag   = G.u_tag;

        [it0, iu0] = find_first_valid_(chi_map);
        T_slider.Limits = [1, max(1,numel(T_list))];
        U_slider.Limits = [1, max(1,numel(U_list))];
        T_slider.Value  = it0;
        U_slider.Value  = iu0;

        psi1_opt = 0; psi2_opt = 0;
        refresh_labels_();
        update_plot_only_();
    end

    function render_empty_(msg)
        cla(ax2d); cla(axPsi); cla(axX);
        text(ax2d, 0.5, 0.5, msg, 'Units','normalized', 'HorizontalAlignment','center', 'FontSize', 14);
        coeff_box.Value = {msg};
    end
end

% =====================================================================
% Embedded loader (same minimal idea as button version, plus mu_folder_map)
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
    if ukind=="mu", muF(k)=uval; end
end

if u_tag=="", u_tag="doping"; end

mask = isfinite(Tvals) & isfinite(Uvals) & isfinite(chiV) & isfinite(dopH);
Tvals=Tvals(mask); Uvals=Uvals(mask); chiV=chiV(mask); dopH=dopH(mask); muF=muF(mask);

T_list = sort(unique(Tvals));
U_list = sort(unique(Uvals));
NT = numel(T_list); NU = numel(U_list);

chi_sum=zeros(NT,NU); chi_cnt=zeros(NT,NU);
dop_sum=zeros(NT,NU); dop_cnt=zeros(NT,NU);
mu_sum=zeros(NT,NU);  mu_cnt=zeros(NT,NU);

for i=1:numel(Tvals)
    iT=find(T_list==Tvals(i),1);
    iU=find(U_list==Uvals(i),1);
    if isempty(iT)||isempty(iU), continue; end
    chi_sum(iT,iU)=chi_sum(iT,iU)+chiV(i); chi_cnt(iT,iU)=chi_cnt(iT,iU)+1;
    dop_sum(iT,iU)=dop_sum(iT,iU)+dopH(i); dop_cnt(iT,iU)=dop_cnt(iT,iU)+1;
    if isfinite(muF(i))
        mu_sum(iT,iU)=mu_sum(iT,iU)+muF(i); mu_cnt(iT,iU)=mu_cnt(iT,iU)+1;
    end
end

chi_map=nan(NT,NU); m=chi_cnt>0; chi_map(m)=chi_sum(m)./chi_cnt(m);
doping_map=nan(NT,NU); md=dop_cnt>0; doping_map(md)=dop_sum(md)./dop_cnt(md);

mu_folder_map=nan(NT,NU); mm=mu_cnt>0; mu_folder_map(mm)=mu_sum(mm)./mu_cnt(mm);

G = struct('T_list',T_list,'U_list',U_list,'chi_map',chi_map,'doping_map',doping_map, ...
           'mu_folder_map',mu_folder_map,'u_tag',u_tag);
end

function [Tfolder, Ufolder, ukind, uval] = parse_TU_from_path_(fpath)
Tfolder=""; Ufolder=""; ukind=""; uval=NaN; %#ok<NASGU>
p = replace(string(fpath), "\", "/");
tokT = regexp(p, '/(T[^/]+)/', 'tokens', 'once'); %#ok<NASGU>
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
fid=fopen(fpath,'r'); if fid<0, return; end
c=onCleanup(@() fclose(fid)); %#ok<NASGU>
num='([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';
for t=1:3000
    ln=fgetl(fid); if ~ischar(ln), break; end
    s=strtrim(ln); if isempty(s), continue; end
    if s(1)~='#', break; end
    if ~isfinite(H.T_K)
        tok=regexp(s,['^\#\s*T\s*=\s*' num],'tokens','once');
        if isempty(tok), tok=regexp(s,['^\#\s*T_K\s*=\s*' num],'tokens','once'); end
        if ~isempty(tok), H.T_K=str2double(tok{1}); end
    end
    if ~isfinite(H.doping)
        tok=regexp(s,['^\#\s*doping\s*=\s*' num],'tokens','once');
        if ~isempty(tok), H.doping=str2double(tok{1}); end
    end
    if isfinite(H.T_K) && isfinite(H.doping), break; end
end
H.ok = isfinite(H.T_K) && isfinite(H.doping);
end

function M = read_numeric_skiphash_(fpath)
fid=fopen(fpath,'r'); if fid<0, M=[]; return; end
c=onCleanup(@() fclose(fid)); %#ok<NASGU>
rows={};
while true
    ln=fgetl(fid); if ~ischar(ln), break; end
    if isempty(ln), continue; end
    if ~isempty(regexp(ln,'^\s*#','once')), continue; end
    v=sscanf(ln,'%f').';
    if isempty(v), continue; end
    rows{end+1,1}=v; %#ok<AGROW>
end
if isempty(rows), M=[]; return; end
ncol=max(cellfun(@numel,rows));
M=nan(numel(rows),ncol);
for i=1:numel(rows)
    v=rows{i}; M(i,1:numel(v))=v;
end
lastFinite=find(any(isfinite(M),1),1,'last');
if ~isempty(lastFinite), M=M(:,1:lastFinite); end
end

function C = detect_cols_(M)
ncol=size(M,2);
C=struct('iq',1,'jq',2,'Re',5);
if ncol>=8, C.iq=2; C.jq=3; C.Re=6; end
end

function [it, iu] = find_first_valid_(M)
[it, iu] = find(isfinite(M), 1, 'first');
if isempty(it), it=1; iu=1; end
end

function i = clamp_(i,N), i=max(1,min(N,i)); end

function iu2 = nearest_valid_in_row_(M, it, iu0)
row=M(it,:);
good=find(isfinite(row));
if isempty(good), iu2=iu0; return; end
[~,k]=min(abs(good-iu0));
iu2=good(k);
end