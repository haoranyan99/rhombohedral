function out_png = plot_fermi_E_and_f_scan_T_mu_UI2()
% plot_fermi_E_and_f_scan_T_mu_UI3
% Browse debug fermiPatch files with THREE sliders:
%   Slider-1: Temperature T (from header; fallback to /Txx/ in path if header missing)
%   Slider-2: mu-folder value (from path token /mu0.985000/; fallback to header EF if missing)
%   Slider-3: file idx under that (T,mu) (keep original scan order; no sorting within group)
%
% Root selected: e.g. data/fermi_sk_mu_200/D0.067
% Scan pattern: root/T*/**/fermiPatch*.txt
%
% Header format (example):
%   # mode = mu
%   # T = 0
%   # EF = 0.945
%   # doping = 0
%   # filling = 0.5
%   # kmesh = BZ Nk=100
%   # b1 = (...)
%   # b2 = (...)
%
% Data columns (debug):
%   idx iq jq kx ky occ_k_avg  E0..E{dim-1}  f0..f{dim-1}
%
% Output:
%   out_png = last saved png path

out_png = "";

% -----------------------------
% choose root directory
% -----------------------------
start_dir = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\fermi_sk_mu_200\D0.067\";
root_dir = uigetdir(start_dir, "Select root directory like data/fermiXXX/D0.067");
if isequal(root_dir,0)
    error("No directory selected.");
end
root_dir = string(root_dir);

% -----------------------------
% scan files
% -----------------------------
fprintf("[scan] root_dir = %s\n", root_dir);
tscan = tic;
files = dir(fullfile(root_dir, "T*", "**", "fermiPatch*.txt"));
fprintf("[scan] dir(T*/**/fermiPatch*.txt) -> %d files in %.3f s\n", numel(files), toc(tscan));
if isempty(files)
    error("No fermiPatch*.txt found under: %s", root_dir);
end

DB = scan_fermi_tree_by_T_mu_keep_order(files);
if isempty(DB.Tvals)
    error("No valid files: need T in header or /Txx/ in path.");
end

% initial indices
it = 1;
imu = 1;
id  = 1;

% -----------------------------
% typography
% -----------------------------
FS = 14;
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
set(groot, "defaultTextInterpreter", "latex");

% -----------------------------
% build figure
% -----------------------------
fig = figure("Color","w", "Units","pixels", "Position",[80 80 1440 760], ...
    "Name","fermiPatch En/fn browser (T + mu + file idx)");

hInfo = uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.05 0.965 0.93 0.030], ...
    "String","", "BackgroundColor","w", ...
    "FontSize", 11, "HorizontalAlignment","left");

tl = tiledlayout(fig, 1, 2, "TileSpacing","compact", "Padding","compact");
tl.Position = [0.06 0.20 0.92 0.73];
ax1 = nexttile(tl, 1);
ax2 = nexttile(tl, 2);

% -----------------------------
% sliders
% -----------------------------
% T slider
uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.06 0.145 0.08 0.03], ...
    "String","T idx", "BackgroundColor","w", "FontSize", 11, ...
    "HorizontalAlignment","left");

hT = uicontrol(fig, "Style","slider", "Units","normalized", ...
    "Position",[0.13 0.148 0.56 0.025], ...
    "Min",1, "Max",numel(DB.Tvals), "Value",it, ...
    "SliderStep", slider_steps(numel(DB.Tvals)), ...
    "Callback", @onTChanged);

tT = uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.70 0.145 0.10 0.03], ...
    "String","", "BackgroundColor","w", "FontSize", 11, ...
    "HorizontalAlignment","left");

% mu slider
uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.06 0.105 0.08 0.03], ...
    "String","mu idx", "BackgroundColor","w", "FontSize", 11, ...
    "HorizontalAlignment","left");

hMu = uicontrol(fig, "Style","slider", "Units","normalized", ...
    "Position",[0.13 0.108 0.56 0.025], ...
    "Min",1, "Max",max(1,numel(DB.byT{it}.muvals)), "Value",imu, ...
    "SliderStep", slider_steps(max(1,numel(DB.byT{it}.muvals))), ...
    "Callback", @onMuChanged);

tMu = uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.70 0.105 0.20 0.03], ...
    "String","", "BackgroundColor","w", "FontSize", 11, ...
    "HorizontalAlignment","left");

% file idx slider
uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.06 0.065 0.18 0.03], ...
    "String","file idx (under T,mu)", "BackgroundColor","w", "FontSize", 11, ...
    "HorizontalAlignment","left");

nfile0 = numel(DB.byT{it}.byMu{imu}.files);
hD = uicontrol(fig, "Style","slider", "Units","normalized", ...
    "Position",[0.24 0.068 0.45 0.025], ...
    "Min",1, "Max",max(1,nfile0), "Value",id, ...
    "SliderStep", slider_steps(max(1,nfile0)), ...
    "Callback", @onFileChanged);

% band selector
uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.78 0.145 0.08 0.03], ...
    "String","band n", "BackgroundColor","w", "FontSize", 11, ...
    "HorizontalAlignment","left");

hBand = uicontrol(fig, "Style","edit", "Units","normalized", ...
    "Position",[0.86 0.148 0.05 0.030], ...
    "String","0", "FontSize", 12, ...
    "Callback", @onBandChanged);

uicontrol(fig, "Style","pushbutton", "Units","normalized", ...
    "Position",[0.92 0.144 0.06 0.038], ...
    "String","Apply", "FontSize", 12, ...
    "Callback", @onBandChanged);

uicontrol(fig, "Style","pushbutton", "Units","normalized", ...
    "Position",[0.78 0.060 0.20 0.050], ...
    "String","Save PNG", "FontSize", 12, ...
    "Callback", @onSave);

% -----------------------------
% state
% -----------------------------
state = struct();
state.root_dir = root_dir;
state.DB = DB;
state.it = it;
state.imu = imu;
state.id = id;
state.n  = int32(0);
state.in_path = "";
state.meta = [];
state.D = [];

% datacursor
dcm = datacursormode(fig);
set(dcm, "Enable","on", "DisplayStyle","window", "SnapToDataVertex","off");
set(dcm, "UpdateFcn", @data_tip_update);

% initial load & plot
reload_current_file();
update_plot();
sync_label_texts();

% =============================
% callbacks
% =============================
function onTChanged(~,~)
    state.it = clamp_int(round(get(hT,"Value")), 1, numel(state.DB.Tvals));
    set(hT, "Value", state.it);

    % reset mu slider range for this T
    nmu = numel(state.DB.byT{state.it}.muvals);
    if nmu < 1, nmu = 1; end
    state.imu = clamp_int(state.imu, 1, nmu);
    set(hMu, "Min",1, "Max",nmu, "Value",state.imu, ...
        "SliderStep", slider_steps(nmu));

    % reset file slider for this (T,mu)
    nx = numel(state.DB.byT{state.it}.byMu{state.imu}.files);
    if nx < 1, nx = 1; end
    state.id = clamp_int(state.id, 1, nx);
    set(hD, "Min",1, "Max",nx, "Value",state.id, ...
        "SliderStep", slider_steps(nx));

    reload_current_file();
    update_plot();
    sync_label_texts();
end

function onMuChanged(~,~)
    nmu = numel(state.DB.byT{state.it}.muvals);
    state.imu = clamp_int(round(get(hMu,"Value")), 1, max(1,nmu));
    set(hMu, "Value", state.imu);

    nx = numel(state.DB.byT{state.it}.byMu{state.imu}.files);
    if nx < 1, nx = 1; end
    state.id = clamp_int(1, 1, nx); % move to first file when changing mu
    set(hD, "Min",1, "Max",nx, "Value",state.id, ...
        "SliderStep", slider_steps(nx));

    reload_current_file();
    update_plot();
    sync_label_texts();
end

function onFileChanged(~,~)
    nx = numel(state.DB.byT{state.it}.byMu{state.imu}.files);
    state.id = clamp_int(round(get(hD,"Value")), 1, max(1,nx));
    set(hD, "Value", state.id);

    reload_current_file();
    update_plot();
    sync_label_texts();
end

function onBandChanged(~,~)
    v = str2double(get(hBand,"String"));
    if ~(isfinite(v) && abs(v-round(v))<1e-12)
        set(hBand,"String", string(state.n));
        return;
    end
    state.n = int32(round(v));
    if ~isempty(state.D) && isfield(state.D,"dim") && state.D.dim > 0
        state.n = int32(clamp_int(double(state.n), 0, state.D.dim-1));
        set(hBand,"String", string(state.n));
    end
    update_plot();
end

function onSave(~,~)
    if strlength(state.in_path) == 0, return; end
    [fdir, base, ~] = fileparts(state.in_path);
    out_png = fullfile(string(fdir), base + "_En_fn_UI_n" + string(state.n) + ".png");
    exportgraphics(fig, out_png, "Resolution", 300);
    fprintf("Saved:\n  %s\n", out_png);
end

% =============================
% core ops
% =============================
function reload_current_file()
    entryT = state.DB.byT{state.it};
    entryM = entryT.byMu{state.imu};

    nx = numel(entryM.files);
    if nx < 1
        state.in_path = "";
        state.meta = [];
        state.D = [];
        set(hInfo,"String", sprintf("No files under T=%gK, mu=%g", entryT.T, entryM.mu));
        return;
    end

    state.id = clamp_int(state.id, 1, nx);
    state.in_path = entryM.files(state.id);

    meta = parse_header_meta_fermi_debug_single_chistyle(state.in_path); % header-only
    D = load_and_process_fermi_debug_single(state.in_path, meta);

    state.meta = meta;
    state.D = D;

    if D.dim > 0
        state.n = int32(clamp_int(double(state.n), 0, D.dim-1));
        set(hBand,"String", string(state.n));
    else
        state.n = int32(0);
        set(hBand,"String", "0");
    end
end

function sync_label_texts()
    Tnow = state.DB.Tvals(state.it);
    set(tT, "String", sprintf("T=%.6g", Tnow));

    muvals = state.DB.byT{state.it}.muvals;
    if isempty(muvals)
        set(tMu, "String", "mu=(none)");
    else
        munow = muvals(state.imu);
        set(tMu, "String", sprintf("\\mu(folder)=%.8g", munow));
    end
end

function update_plot()
    cla(ax1); cla(ax2);

    if isempty(state.D) || ~state.D.has_debug_bands
        text(ax1, 0.5, 0.5, "No debug bands (E/f) in file", "Units","normalized", ...
            "HorizontalAlignment","center", "FontSize", FS);
        text(ax2, 0.5, 0.5, "No debug bands (E/f) in file", "Units","normalized", ...
            "HorizontalAlignment","center", "FontSize", FS);
        set(hInfo,"String", "No debug bands (E/f).");
        return;
    end

    n = double(state.n) + 1;
    EG = state.D.EG(:,:,n);
    fG = state.D.fG(:,:,n);

    maskE = isfinite(EG) & isfinite(state.D.kx0G) & isfinite(state.D.ky0G);
    maskF = isfinite(fG) & isfinite(state.D.kx0G) & isfinite(state.D.ky0G);

    dop0  = state.D.meta.doping;
    fill0 = state.D.meta.filling;
    EF0   = state.D.meta.EF_eV;
    T0    = state.D.meta.T_K;
    mode0 = "";
    if isfield(state.D.meta,"mode"), mode0 = string(state.D.meta.mode); end

    Emin = NaN; Emax = NaN;
    fmin = NaN; fmax = NaN;
    if any(maskE,'all'), Emin = min(EG(maskE)); Emax = max(EG(maskE)); end
    if any(maskF,'all'), fmin = min(fG(maskF)); fmax = max(fG(maskF)); end

    % show also which (T,mu,file)
    entryT = state.DB.byT{state.it};
    entryM = entryT.byMu{state.imu};
    set(hInfo, "String", sprintf( ...
        "T=%.6g K | mu(folder)=%.10g | file idx=%d/%d | mode=%s | EF(mu)=%.10g eV | doping=%.10g | filling=%.16g | Nk=%g | band n=%d | En[min,max]=[%.6g,%.6g] | fn[min,max]=[%.6g,%.6g]", ...
        T0, entryM.mu, state.id, numel(entryM.files), char(mode0), EF0, dop0, fill0, state.D.meta.Nk, state.n, Emin, Emax, fmin, fmax));

    % LEFT: En
    Z1 = EG; Z1(~maskE) = NaN;
    hS1 = surf(ax1, state.D.kx0G, state.D.ky0G, Z1, "EdgeColor","none");
    view(ax1, 2); axis(ax1,"equal"); axis(ax1,"tight");
    xlabel(ax1, '$(k_x-k_x^0)$', "FontSize", FS);
    ylabel(ax1, '$(k_y-k_y^0)$', "FontSize", FS);
    title(ax1, sprintf('$E_{%d}(\\mathbf{k})$ (eV)', state.n), "FontSize", FS, "FontWeight","normal");
    set(ax1, "FontSize", FS, "LineWidth", 1.0, "TickDir","out", "Box","on");
    colormap(ax1, jet(256));
    cb1 = colorbar(ax1, 'northoutside'); cb1.Box = 'on';
    hold(ax1,"on"); draw_BZ_rhombus_boundary(ax1, state.D); hold(ax1,"off");
    setappdata(hS1, "panel_tag", "E");

    % RIGHT: fn
    Z2 = fG; Z2(~maskF) = NaN;
    hS2 = surf(ax2, state.D.kx0G, state.D.ky0G, Z2, "EdgeColor","none");
    view(ax2, 2); axis(ax2,"equal"); axis(ax2,"tight");
    xlabel(ax2, '$(k_x-k_x^0)$', "FontSize", FS);
    ylabel(ax2, '$(k_y-k_y^0)$', "FontSize", FS);
    title(ax2, sprintf('$f_{%d}(\\mathbf{k})$', state.n), "FontSize", FS, "FontWeight","normal");
    set(ax2, "FontSize", FS, "LineWidth", 1.0, "TickDir","out", "Box","on");
    colormap(ax2, flipud(hot(256))); caxis(ax2, [0,1]);
    cb2 = colorbar(ax2, 'northoutside'); cb2.Box = 'on';
    hold(ax2,"on"); draw_BZ_rhombus_boundary(ax2, state.D); hold(ax2,"off");
    setappdata(hS2, "panel_tag", "f");
end

% =============================
% datatip
% =============================
function txt = data_tip_update(~, event)
    if isempty(state.D) || ~state.D.has_debug_bands
        txt = {"(no data)"}; return;
    end

    pos = event.Position; xq = pos(1); yq = pos(2);
    X = state.D.kx0G; Y = state.D.ky0G;

    good = isfinite(X) & isfinite(Y);
    dx = X - xq; dy = Y - yq;
    dist2 = dx.^2 + dy.^2; dist2(~good) = inf;
    [~, lin] = min(dist2(:));
    [rr, cc] = ind2sub(size(dist2), lin);

    idx0 = state.D.idxG(rr,cc);
    iq0  = state.D.iqG(rr,cc);
    jq0  = state.D.jqG(rr,cc);

    kx0 = X(rr,cc); ky0 = Y(rr,cc);
    n = double(state.n) + 1;
    En = state.D.EG(rr,cc,n);
    fn = state.D.fG(rr,cc,n);

    tag = "";
    try, tag = string(getappdata(event.Target, "panel_tag")); catch, end

    txt = {
        sprintf("idx = %g", idx0)
        sprintf("(iq,jq) = (%g, %g)", iq0, jq0)
        sprintf("(kx,ky) = (%.8g, %.8g)", kx0, ky0)
        sprintf("band n = %d", state.n)
        sprintf("E_n = %.10g eV", En)
        sprintf("f_n = %.10g", fn)
    };
    if strlength(tag) > 0, txt{end+1} = sprintf("panel = %s", tag); end
end

end % ===== end main =====


% =========================================================
% DB scan: group by T, then by mu(folder), KEEP original order
% =========================================================
function DB = scan_fermi_tree_by_T_mu_keep_order(files)
rec_T  = [];
rec_mu = [];
rec_p  = strings(0,1);

nfile = numel(files);
for i = 1:nfile
    if mod(i,200)==1
        fprintf("[scan] header %d/%d\n", i, nfile);
    end

    fpath = string(fullfile(files(i).folder, files(i).name));
    fp = strrep(fpath, "\", "/");

    % --- T: header first, else /Txx/
    m2 = parse_header_T_doping_mu_like_chi(fpath);
    Tval = m2.T;
    if ~isfinite(Tval)
        tokT = regexp(fp, "/T([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)/", "tokens", "once");
        if ~isempty(tokT), Tval = str2double(tokT{1}); end
    end
    if ~isfinite(Tval), continue; end

    % --- mu: prefer folder token /mu0.985000/ ; fallback to header mu(EF)
    muval = parse_mu_from_path(fp);
    if ~isfinite(muval)
        muval = m2.mu; % EF from header
    end
    if ~isfinite(muval)
        muval = NaN; % still allow grouping as NaN? we skip to avoid mess
        continue;
    end

    rec_T(end+1,1)  = Tval; %#ok<AGROW>
    rec_mu(end+1,1) = muval; %#ok<AGROW>
    rec_p(end+1,1)  = fpath; %#ok<AGROW>
end

if isempty(rec_T)
    DB = struct("Tvals",[],"byT",{{}});
    return;
end

Ttol  = 1e-9;
mutol = 1e-12;

Tvals = uniquetol(rec_T, Ttol);
Tvals = sort(Tvals);

byT = cell(numel(Tvals),1);
for it = 1:numel(Tvals)
    Tv = Tvals(it);
    maskT = abs(rec_T - Tv) < Ttol;

    mu_here = rec_mu(maskT);
    p_here  = rec_p(maskT);

    muvals = uniquetol(mu_here, mutol);
    muvals = sort(muvals);

    byMu = cell(numel(muvals),1);
    for j = 1:numel(muvals)
        muv = muvals(j);
        maskM = abs(mu_here - muv) < mutol;

        % keep original order
        ps = p_here(maskM);

        byMu{j} = struct();
        byMu{j}.mu    = muv;
        byMu{j}.files = ps(:);
        byMu{j}.nfile = numel(ps);
    end

    byT{it} = struct();
    byT{it}.T     = Tv;
    byT{it}.muvals= muvals(:);
    byT{it}.byMu  = byMu;
end

DB = struct();
DB.Tvals = Tvals(:);
DB.byT   = byT;
end

function mu = parse_mu_from_path(fp)
% parse folder token ".../mu0.985000/..."
mu = NaN;
tok = regexp(fp, "/mu([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)(?:/|$)", "tokens", "once");
if ~isempty(tok)
    mu = str2double(tok{1});
end
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
% Full meta parse for plotting (still header-only, chi-style)
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
    if ~isfinite(meta.dk_frac)
        tokDk = regexp(s, "(?:^|\s)dk_frac\s*=\s*"+num, "tokens", "once");
        if ~isempty(tokDk), meta.dk_frac = str2double(tokDk{1}); end
    end

    if ~isfinite(meta.T_K)
        tok = regexp(s, "(?:^|\s)T_K\s*=\s*"+num, "tokens", "once");
        if ~isempty(tok), meta.T_K = str2double(tok{1}); end
    end
    if ~isfinite(meta.EF_eV)
        tok = regexp(s, "(?:^|\s)EF_eV\s*=\s*"+num, "tokens", "once");
        if ~isempty(tok), meta.EF_eV = str2double(tok{1}); end
    end
end
end

% =========================================================
% loader: build (2Nk+1)x(2Nk+1) grid and E/f cubes
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
    fprintf("[En/fn UI] meta.Nk missing -> inferred Nk=%d\n", meta.Nk);
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
    EG = NaN(Nk_side, Nk_side, dim);
    fG = NaN(Nk_side, Nk_side, dim);
else
    EG = [];
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
            Erow = raw(r, 7 : 6+dim);
            frow = raw(r, 7+dim : 6+2*dim);
            EG(rr,cc,:) = reshape(Erow, 1,1,dim);
            fG(rr,cc,:) = reshape(frow, 1,1,dim);
        end
    end
end

% origin shift by (iq,jq)=(0,0) if exists, else first available
r0 = 0 + Nk_rad + 1;
c0 = 0 + Nk_rad + 1;
if ~(isfinite(kxG(r0,c0)) && isfinite(kyG(r0,c0)))
    [rr,cc] = find(has, 1, "first");
    r0 = rr; c0 = cc;
    fprintf("[En/fn UI] WARNING: (iq,jq)=(0,0) missing; using first available cell as origin.\n");
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
D.EG  = EG;
D.fG  = fG;
end

function draw_BZ_rhombus_boundary(ax, D)
if ~isfield(D, "meta") || ~isfield(D.meta, "b1") || ~isfield(D.meta, "b2")
    return;
end
b1 = D.meta.b1;
b2 = D.meta.b2;
if any(~isfinite([b1; b2]))
    return;
end

k0 = [0;0];
if isfield(D, "k0") && all(isfinite(D.k0))
    k0 = D.k0;
end

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