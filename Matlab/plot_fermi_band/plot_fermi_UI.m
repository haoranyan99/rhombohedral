function out = plot_fermi_UI()
% plot_fermi_UI_copyfig
% Minimal UI:
%   - Only plot f_n(k) heatmap (NO energy panel)
%   - Two sliders: T idx, X idx (X = mu(folder preferred) OR doping(header))
%   - Manual ranges at top (T_range, mu_range/dop_range)
%   - band_n fixed by user setting (not adjustable)
%   - Button: "Copy Figure" duplicates current view to a NEW figure window
%
% Root: .../fermi_sk_mu_200/D0.067
% Scan: root/T*/**/fermiPatch*.txt

out = struct();
out.png_main = "";
out.png_last_copy = "";

% ============================================================
% USER SETTINGS (edit here)
% ============================================================
start_dir = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\fermi_sk_mu_200\D0.084\";

% X_mode = "mu" or "doping"
X_mode = "mu";

% ---- manual ranges ( [] => no filter ) ----
T_range  = [4.0, 4.0];        % e.g. [3.5,4.5] or [] or [4,4]
mu_range = [1.00, 1e5];     % ONLY used when X_mode="mu"
dop_range= [0.9,  1.1];       % ONLY used when X_mode="doping"

band_n = 6;                   % fixed band index (0-based)
pick_policy = "newest";       % "newest" or "oldest" per (T,X)

FS = 14;

% ============================================================
% choose root directory
% ============================================================
root_dir = uigetdir(start_dir, "Select root directory like data/fermiXXX/D0.067");
if isequal(root_dir,0), error("No directory selected."); end
root_dir = string(root_dir);

% ============================================================
% scan files
% ============================================================
fprintf("[scan] root_dir = %s\n", root_dir);
tscan = tic;
files = dir(fullfile(root_dir, "T*", "**", "fermiPatch*.txt"));
fprintf("[scan] dir(T*/**/fermiPatch*.txt) -> %d files in %.3f s\n", numel(files), toc(tscan));
if isempty(files), error("No fermiPatch*.txt found under: %s", root_dir); end

DB = scan_fermi_tree_by_T_X_keep_pick(files, X_mode, T_range, mu_range, dop_range, pick_policy);
if isempty(DB.Tvals)
    error("No valid files after filtering. Check T_range / mu_range / dop_range / headers.");
end

% initial indices
it = 1;
ix = 1;

% typography
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
set(groot, "defaultTextInterpreter", "latex");

% ============================================================
% build main figure (single panel)
% ============================================================
fig = figure("Color","w", "Units","pixels", "Position",[80 80 1200 760], ...
    "Name","fermiPatch f_n(k) browser (Copy Figure)");

hInfo = uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.05 0.955 0.93 0.035], ...
    "String","", "BackgroundColor","w", ...
    "FontSize", 11, "HorizontalAlignment","left");

axF = axes(fig, "Position",[0.08 0.22 0.86 0.70]);

% ---- sliders + numeric labels ----
uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.08 0.145 0.05 0.03], ...
    "String","T", "BackgroundColor","w", "FontSize", 11, ...
    "HorizontalAlignment","left");

hT = uicontrol(fig, "Style","slider", "Units","normalized", ...
    "Position",[0.13 0.148 0.60 0.025], ...
    "Min",1, "Max",numel(DB.Tvals), "Value",it, ...
    "SliderStep", slider_steps(numel(DB.Tvals)), ...
    "Callback", @onTChanged);

tT = uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.74 0.145 0.22 0.03], ...
    "String","", "BackgroundColor","w", "FontSize", 11, ...
    "HorizontalAlignment","left");

uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.08 0.105 0.05 0.03], ...
    "String", char(X_mode), "BackgroundColor","w", "FontSize", 11, ...
    "HorizontalAlignment","left");

hX = uicontrol(fig, "Style","slider", "Units","normalized", ...
    "Position",[0.13 0.108 0.60 0.025], ...
    "Min",1, "Max",max(1,numel(DB.byT{it}.Xvals)), "Value",ix, ...
    "SliderStep", slider_steps(max(1,numel(DB.byT{it}.Xvals))), ...
    "Callback", @onXChanged);

tX = uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.74 0.105 0.22 0.03], ...
    "String","", "BackgroundColor","w", "FontSize", 11, ...
    "HorizontalAlignment","left");

% ---- buttons ----
uicontrol(fig, "Style","pushbutton", "Units","normalized", ...
    "Position",[0.08 0.045 0.16 0.055], ...
    "String","Copy Figure", "FontSize", 12, ...
    "Callback", @onCopyFigure);

uicontrol(fig, "Style","pushbutton", "Units","normalized", ...
    "Position",[0.26 0.045 0.16 0.055], ...
    "String","Save PNG", "FontSize", 12, ...
    "Callback", @onSaveMain);

% ============================================================
% state
% ============================================================
state = struct();
state.root_dir = root_dir;
state.DB = DB;
state.it = it;
state.ix = ix;
state.band_n = band_n;
state.in_path = "";
state.meta = [];
state.D = [];

% datacursor
dcm = datacursormode(fig);
set(dcm, "Enable","on", "DisplayStyle","window", "SnapToDataVertex","off");
set(dcm, "UpdateFcn", @data_tip_update);

% initial
reload_current_file();
sync_labels();
update_plot();

% ============================================================
% callbacks
% ============================================================
function onTChanged(~,~)
    state.it = clamp_int(round(get(hT,"Value")), 1, numel(state.DB.Tvals));
    set(hT, "Value", state.it);

    % reset X slider for this T
    nX = numel(state.DB.byT{state.it}.Xvals);
    if nX < 1, nX = 1; end
    state.ix = clamp_int(state.ix, 1, nX);
    set(hX, "Min",1, "Max",nX, "Value",state.ix, ...
        "SliderStep", slider_steps(nX));

    reload_current_file();
    sync_labels();
    update_plot();
end

function onXChanged(~,~)
    nX = numel(state.DB.byT{state.it}.Xvals);
    state.ix = clamp_int(round(get(hX,"Value")), 1, max(1,nX));
    set(hX, "Value", state.ix);

    reload_current_file();
    sync_labels();
    update_plot();
end

function onCopyFigure(~,~)
    if isempty(state.D) || ~state.D.has_debug_bands
        fprintf("[Copy] No valid debug bands to copy.\n");
        return;
    end

    % make a new figure and replot the same content
    fig2 = figure("Color","w","Units","pixels","Position",[140 140 900 720], ...
        "Name", make_copy_title());
    ax2 = axes(fig2, "Position",[0.10 0.12 0.83 0.80]);

    % plot f on new axes
    plot_f_on_axes(ax2);

    drawnow;
end

function s = make_copy_title()
    [Tnow, Xnow] = current_TX();
    if X_mode=="mu"
        s = sprintf("Copy: T=%.6g, mu=%.10g, band_n=%d", Tnow, Xnow, state.band_n);
    else
        s = sprintf("Copy: T=%.6g, doping=%.10g, band_n=%d", Tnow, Xnow, state.band_n);
    end
end

function onSaveMain(~,~)
    if strlength(state.in_path)==0, return; end
    [fdir, base, ~] = fileparts(state.in_path);
    tag = sprintf("_Tidx%d_Xidx%d_n%d", state.it, state.ix, state.band_n);
    out.png_main = fullfile(string(fdir), base + tag + "_f.png");
    exportgraphics(fig, out.png_main, "Resolution", 300);
    fprintf("Saved main:\n  %s\n", out.png_main);
end

% ============================================================
% core
% ============================================================
function [Tnow, Xnow] = current_TX()
    Tnow = state.DB.Tvals(state.it);
    Xvals = state.DB.byT{state.it}.Xvals;
    if isempty(Xvals), Xnow = NaN;
    else, Xnow = Xvals(state.ix);
    end
end

function reload_current_file()
    entryT = state.DB.byT{state.it};
    if isempty(entryT.byX) || state.ix > numel(entryT.byX)
        state.in_path = "";
        state.meta = [];
        state.D = [];
        return;
    end

    entryX = entryT.byX{state.ix};
    if strlength(entryX.pick) == 0
        state.in_path = "";
        state.meta = [];
        state.D = [];
        return;
    end

    state.in_path = entryX.pick;

    meta = parse_header_meta_fermi_debug_single_chistyle(state.in_path);
    D = load_and_process_fermi_debug_single(state.in_path, meta);

    state.meta = meta;
    state.D = D;
end

function sync_labels()
    Tnow = state.DB.Tvals(state.it);
    set(tT, "String", sprintf("T=%.6g", Tnow));

    Xvals = state.DB.byT{state.it}.Xvals;
    if isempty(Xvals)
        set(tX, "String", sprintf("%s=(none)", char(X_mode)));
    else
        Xnow = Xvals(state.ix);
        if X_mode == "mu"
            set(tX, "String", sprintf("\\mu=%.10g", Xnow));
        else
            set(tX, "String", sprintf("doping=%.10g", Xnow));
        end
    end
end

function update_plot()
    cla(axF);

    if isempty(state.D) || ~state.D.has_debug_bands
        text(axF, 0.5, 0.5, "No debug bands (f) in file", "Units","normalized", ...
            "HorizontalAlignment","center", "FontSize", FS);
        axis(axF,"off");
        set(hInfo,"String","No debug bands (f).");
        return;
    end

    % fixed band index
    n = double(state.band_n) + 1;
    if n < 1 || n > state.D.dim
        text(axF, 0.5, 0.5, sprintf("band_n=%d out of range [0,%d]", state.band_n, state.D.dim-1), ...
            "Units","normalized", "HorizontalAlignment","center", "FontSize", FS);
        axis(axF,"off");
        set(hInfo,"String","Band index out of range.");
        return;
    end

    plot_f_on_axes(axF);
end

function plot_f_on_axes(ax)
    n = double(state.band_n) + 1;

    fG = state.D.fG(:,:,n);
    maskF = isfinite(fG) & isfinite(state.D.kx0G) & isfinite(state.D.ky0G);

    Z = fG; Z(~maskF) = NaN;

    surf(ax, state.D.kx0G, state.D.ky0G, Z, "EdgeColor","none");
    view(ax, 2); axis(ax,"equal"); axis(ax,"tight");
    xlabel(ax, '$(k_x-k_x^0)$', "FontSize", FS);
    ylabel(ax, '$(k_y-k_y^0)$', "FontSize", FS);

    % info line
    dop0  = state.D.meta.doping;
    fill0 = state.D.meta.filling;
    EF0   = state.D.meta.EF_eV;
    T0    = state.D.meta.T_K;
    mode0 = "";
    if isfield(state.D.meta,"mode"), mode0 = string(state.D.meta.mode); end
% 
%     title(ax, sprintf('$f_{%d}(\\mathbf{k})$ | T=%.6g | mode=%s | EF=%.6g | dop=%.6g | fill=%.6g', ...
%         state.band_n, T0, char(mode0), EF0, dop0, fill0), ...
%         "FontSize", FS, "FontWeight","normal");

    set(ax, "FontSize", FS, "LineWidth", 1.0, "TickDir","out", "Box","on");
    colormap(ax, flipud(hot(256))); caxis(ax, [0,1]);
    cb = colorbar(ax, 'northoutside'); cb.Box = 'on';

    hold(ax,"on"); draw_BZ_rhombus_boundary(ax, state.D); hold(ax,"off");

    set(hInfo, "String", sprintf("file=%s", char(state.in_path)));
end

% =============================
% datatip
% =============================
function txt = data_tip_update(~, event)
    if isempty(state.D) || ~state.D.has_debug_bands
        txt = {"(no data)"}; return;
    end
    n = double(state.band_n) + 1;
    if n < 1 || n > state.D.dim
        txt = {sprintf("band_n=%d out of range", state.band_n)}; return;
    end

    pos = event.Position; xq = pos(1); yq = pos(2);
    Xg = state.D.kx0G; Yg = state.D.ky0G;

    good = isfinite(Xg) & isfinite(Yg);
    dx = Xg - xq; dy = Yg - yq;
    dist2 = dx.^2 + dy.^2; dist2(~good) = inf;
    [~, lin] = min(dist2(:));
    [rr, cc] = ind2sub(size(dist2), lin);

    idx0 = state.D.idxG(rr,cc);
    iq0  = state.D.iqG(rr,cc);
    jq0  = state.D.jqG(rr,cc);

    kx0 = Xg(rr,cc); ky0 = Yg(rr,cc);
    fn  = state.D.fG(rr,cc,n);

    txt = {
        sprintf("idx = %g", idx0)
        sprintf("(iq,jq) = (%g, %g)", iq0, jq0)
        sprintf("(kx,ky) = (%.8g, %.8g)", kx0, ky0)
        sprintf("band n = %d", state.band_n)
        sprintf("f_n = %.10g", fn)
    };
end

end % ===== end main =====


% =========================================================
% DB scan: group by T then by X; pick one file per group
% =========================================================
function DB = scan_fermi_tree_by_T_X_keep_pick(files, X_mode, T_range, mu_range, dop_range, pick_policy)

rec_T = [];
rec_X = [];
rec_p = strings(0,1);
rec_dn= [];

nfile = numel(files);
for i = 1:nfile
    if mod(i,200)==1
        fprintf("[scan] header %d/%d\n", i, nfile);
    end

    fpath = string(fullfile(files(i).folder, files(i).name));
    fp = strrep(fpath, "\", "/");

    m2 = parse_header_T_doping_mu_like_chi(fpath);

    % --- T: header first, else /Txx/
    Tval = m2.T;
    if ~isfinite(Tval)
        tokT = regexp(fp, "/T([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)/", "tokens", "once");
        if ~isempty(tokT), Tval = str2double(tokT{1}); end
    end
    if ~isfinite(Tval), continue; end
    if ~pass_range(Tval, T_range), continue; end

    % --- X value + filter
    if X_mode == "mu"
        Xval = parse_mu_from_path(fp);
        if ~isfinite(Xval), Xval = m2.mu; end
        if ~isfinite(Xval), continue; end
        if ~pass_range(Xval, mu_range), continue; end
    else
        Xval = m2.doping;
        if ~isfinite(Xval), continue; end
        if ~pass_range(Xval, dop_range), continue; end
    end

    rec_T(end+1,1)  = Tval; %#ok<AGROW>
    rec_X(end+1,1)  = Xval; %#ok<AGROW>
    rec_p(end+1,1)  = fpath; %#ok<AGROW>
    rec_dn(end+1,1) = files(i).datenum; %#ok<AGROW>
end

if isempty(rec_T)
    DB = struct("Tvals",[],"byT",{{}});
    return;
end

Ttol  = 1e-9;
Xtol  = 1e-12;

Tvals = uniquetol(rec_T, Ttol);
Tvals = sort(Tvals);

byT = cell(numel(Tvals),1);
for it = 1:numel(Tvals)
    Tv = Tvals(it);
    maskT = abs(rec_T - Tv) < Ttol;

    X_here  = rec_X(maskT);
    p_here  = rec_p(maskT);
    dn_here = rec_dn(maskT);

    Xvals = uniquetol(X_here, Xtol);
    Xvals = sort(Xvals);

    byX = cell(numel(Xvals),1);
    for j = 1:numel(Xvals)
        Xv = Xvals(j);
        maskX = abs(X_here - Xv) < Xtol;

        ps = p_here(maskX);
        dns= dn_here(maskX);

        if isempty(ps)
            pick = "";
        else
            if pick_policy == "oldest"
                [~, ii] = min(dns);
            else
                [~, ii] = max(dns);
            end
            pick = ps(ii);
        end

        byX{j} = struct("X",Xv,"pick",pick);
    end

    byT{it} = struct("T",Tv,"Xvals",Xvals(:),"byX",{byX});
end

DB = struct();
DB.Tvals = Tvals(:);
DB.byT   = byT;
end

function tf = pass_range(v, r)
if isempty(r), tf = true; return; end
if ~isfinite(v), tf = false; return; end
lo = min(r); hi = max(r);
tf = (v >= lo) && (v <= hi);
end

function mu = parse_mu_from_path(fp)
mu = NaN;
tok = regexp(fp, "/mu([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)(?:/|$)", "tokens", "once");
if ~isempty(tok), mu = str2double(tok{1}); end
end

function st = slider_steps(N)
if N <= 1, st = [1 1];
else, st = [1/(N-1), min(10/(N-1),1)];
end
end

function x = clamp_int(x, lo, hi)
x = max(lo, min(hi, x));
end

% =========================================================
% chi-style header parsing: only leading '#' lines
% =========================================================
function m = parse_header_T_doping_mu_like_chi(in_path)
m = struct("T",NaN,"doping",NaN,"mu",NaN);

fid = fopen(in_path,'r');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

while true
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    s0 = strtrim(tline);
    if ~startsWith(s0, "#"), break; end
    s = strtrim(erase(s0, "#"));

    if ~isfinite(m.T)
        tok = regexp(s, "(?:^|\s)T(?:_K)?\s*=\s*" + num, "tokens", "once");
        if ~isempty(tok), m.T = str2double(tok{1}); end
    end
    if ~isfinite(m.doping)
        tok = regexp(s, "doping\s*=\s*" + num, "tokens", "once");
        if ~isempty(tok), m.doping = str2double(tok{1}); end
    end
    if ~isfinite(m.mu)
        tok = regexp(s, "(?:^|\s)(?:mu(?:_eV)?|EF(?:_eV)?)\s*=\s*" + num, "tokens", "once");
        if ~isempty(tok), m.mu = str2double(tok{1}); end
    end
end
end

% =========================================================
% Full meta parse for plotting (header-only, chi-style)
% =========================================================
function meta = parse_header_meta_fermi_debug_single_chistyle(in_path)
meta = struct();
meta.Nk      = NaN;
meta.dk_frac = NaN;
meta.filling = NaN;
meta.doping  = NaN;
meta.T_K     = NaN;
meta.Dfield  = NaN;
meta.EF_eV   = NaN;
meta.mode    = "";
meta.b1      = [NaN; NaN];
meta.b2      = [NaN; NaN];

m2 = parse_header_T_doping_mu_like_chi(in_path);
meta.T_K    = m2.T;
meta.doping = m2.doping;
meta.EF_eV  = m2.mu;

fid = fopen(in_path,'r');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

while true
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    s0 = strtrim(tline);
    if ~startsWith(s0,"#"), break; end
    s = strtrim(erase(s0,"#"));

    if strlength(string(meta.mode)) == 0
        tok = regexp(s, "(?:^|\s)mode\s*=\s*([^\r\n#]+)", "tokens", "once");
        if ~isempty(tok)
            meta.mode = string(strtrim(tok{1}));
            meta.mode = strip(meta.mode, '"');
        end
    end

    if any(~isfinite(meta.b1))
        tb1 = regexp(s, "b1\s*=\s*\(\s*"+num+"\s*,\s*"+num+"\s*\)", "tokens", "once");
        if ~isempty(tb1), meta.b1 = [str2double(tb1{1}); str2double(tb1{2})]; end
    end
    if any(~isfinite(meta.b2))
        tb2 = regexp(s, "b2\s*=\s*\(\s*"+num+"\s*,\s*"+num+"\s*\)", "tokens", "once");
        if ~isempty(tb2), meta.b2 = [str2double(tb2{1}); str2double(tb2{2})]; end
    end

    if ~isfinite(meta.filling)
        tok = regexp(s, "(?:^|\s)(?:filling|fill|filling_target)\s*=\s*"+num, "tokens", "once");
        if ~isempty(tok), meta.filling = str2double(tok{1}); end
    end

    if ~isfinite(meta.Dfield)
        tok = regexp(s, "(?:^|\s)(?:Dfield_eV|Dfield)\s*=\s*"+num, "tokens", "once");
        if ~isempty(tok), meta.Dfield = str2double(tok{1}); end
    end

    if ~isfinite(meta.Nk)
        tokNk = regexp(s, "(?:^|\s)Nk\s*=\s*(\d+)", "tokens", "once");
        if ~isempty(tokNk), meta.Nk = str2double(tokNk{1}); end
    end
end
end

% =========================================================
% loader: build grid and E/f cubes
% =========================================================
function D = load_and_process_fermi_debug_single(in_path, meta)
raw = readmatrix(in_path, "FileType","text", "CommentStyle","#");
if isempty(raw) || size(raw,2) < 6
    error("Bad data format: expected >=6 cols in %s", in_path);
end

idx = raw(:,1);
iq  = raw(:,2);
jq  = raw(:,3);
kx  = raw(:,4);
ky  = raw(:,5);

ncol = size(raw,2);
has_debug = (ncol > 6);

dim = 0;
if has_debug
    if mod(ncol-6,2) ~= 0
        error("Debug columns malformed: expected ncol = 6 + 2*dim, got %d", ncol);
    end
    dim = (ncol-6)/2;
    if dim <= 0, has_debug = false; end
end

if isnan(meta.Nk)
    meta.Nk = max([abs(iq(:)); abs(jq(:))]);
    fprintf("[f UI] meta.Nk missing -> inferred Nk=%d\n", meta.Nk);
end
Nk_rad  = meta.Nk;
Nk_side = 2*Nk_rad + 1;

kxG  = NaN(Nk_side, Nk_side);
kyG  = NaN(Nk_side, Nk_side);
idxG = NaN(Nk_side, Nk_side);
iqG  = NaN(Nk_side, Nk_side);
jqG  = NaN(Nk_side, Nk_side);
has  = false(Nk_side, Nk_side);

if has_debug
    EG = NaN(Nk_side, Nk_side, dim); %#ok<NASGU>
    fG = NaN(Nk_side, Nk_side, dim);
else
    fG = [];
end

col = iq + Nk_rad + 1;
row = jq + Nk_rad + 1;
ok = (row>=1 & row<=Nk_side & col>=1 & col<=Nk_side);
idxs = find(ok);

for t = 1:numel(idxs)
    r = idxs(t);
    rr = row(r); cc = col(r);
    if ~has(rr,cc)
        kxG(rr,cc)  = kx(r);
        kyG(rr,cc)  = ky(r);
        idxG(rr,cc) = idx(r);
        iqG(rr,cc)  = iq(r);
        jqG(rr,cc)  = jq(r);
        has(rr,cc)  = true;

        if has_debug
            frow = raw(r, 7+dim : 6+2*dim);
            fG(rr,cc,:) = reshape(frow, 1,1,dim);
        end
    end
end

% origin shift
r0 = 0 + Nk_rad + 1;
c0 = 0 + Nk_rad + 1;
if ~(isfinite(kxG(r0,c0)) && isfinite(kyG(r0,c0)))
    [rr,cc] = find(has, 1, "first");
    r0 = rr; c0 = cc;
end
k0 = [kxG(r0,c0); kyG(r0,c0)];
kx0G = kxG - k0(1);
ky0G = kyG - k0(2);

D = struct();
D.meta = meta;
D.k0 = k0;
D.Nk_rad = Nk_rad;
D.kx0G = kx0G;
D.ky0G = ky0G;
D.idxG = idxG;
D.iqG  = iqG;
D.jqG  = jqG;
D.has_debug_bands = has_debug;
D.dim = dim;
D.fG  = fG;
end

function draw_BZ_rhombus_boundary(ax, D)
if ~isfield(D, "meta") || ~isfield(D.meta, "b1") || ~isfield(D.meta, "b2")
    return;
end
b1 = D.meta.b1; b2 = D.meta.b2;
if any(~isfinite([b1; b2])), return; end

k0 = [0;0];
if isfield(D, "k0") && all(isfinite(D.k0)), k0 = D.k0; end

v1 = (-b1 - b2) / 2;
v2 = ( +b1 - b2) / 2;
v3 = ( +b1 + b2) / 2;
v4 = (-b1 + b2) / 2;

V = [v1(1) v1(2);
     v2(1) v2(2);
     v3(1) v3(2);
     v4(1) v4(2);
     v1(1) v1(2)];

V(:,1) = V(:,1) - k0(1);
V(:,2) = V(:,2) - k0(2);

plot(ax, V(:,1), V(:,2), 'k-', "LineWidth", 1.2, "HandleVisibility","off");
end