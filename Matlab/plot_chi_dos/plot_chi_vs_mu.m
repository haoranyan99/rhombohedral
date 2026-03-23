function out_png = plot_chi_vs_mu()
% plot_chi_vs_mu_only
% Scan chi*.txt recursively under a D folder, parse header T/doping/mu(EF),
% also parse polar from folder name like: polar_meV0.000
% Deduplicate by (T,doping,mu,polar) keeping the OLDEST file.
%
% Build ONE figure:
%   - Re chi(iq,jq) vs mu
% with:
%   - T slider
%   - polar slider
%   - iq / jq edit boxes
%
% File numeric table formats supported:
%   if >=8 cols: idx iq jq qx qy Re Im ...
%   else:        iq jq qx qy Re Im ...

out_png = "";

% -----------------------------
% default iq/jq
% -----------------------------
iq0 = -267;
jq0 = -267;
iq0 = -131;
jq0 = -131;

% -----------------------------
% choose root folder
% -----------------------------
default_root = "/Users/haoranyan/rg_master/data/";
% default_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_200\D0.067";

if ~isfolder(default_root)
    warning("Default folder not found: %s\nFallback to pwd.", default_root);
    default_root = string(pwd);
end

root = uigetdir(default_root, 'Select root folder that CONTAINS chi*.txt (recursive)');
if isequal(root, 0), error('User cancelled.'); end
root = string(root);
fprintf("Root folder:\n  %s\n", root);

% -----------------------------
% resolve where T folders live
% -----------------------------
[chi_run_dir, ~] = resolve_T_root(root);

% -----------------------------
% scan all chi files (recursive)
% -----------------------------
files = dir(fullfile(chi_run_dir, "**", "chi*.txt"));
files = files(~[files.isdir]);
Nraw = numel(files);
if Nraw == 0
    error("No chi*.txt found under: %s", chi_run_dir);
end
fprintf("[scan] found %d raw chi files\n", Nraw);

% -----------------------------
% parse headers and keep OLDEST per (T,doping,mu,polar)
% -----------------------------
meta_map = containers.Map('KeyType','char','ValueType','any');

t0 = tic;
t_last_print = 0;
print_every_sec = 0.5;

n_skip_noT = 0;

for k = 1:Nraw
    p = string(fullfile(files(k).folder, files(k).name));

    m = parse_header_T_doping_mu_mode(p);
    pol = parse_polar_from_path(p);

    if ~isfinite(m.T)
        n_skip_noT = n_skip_noT + 1;
        t_now = toc(t0);
        if (t_now - t_last_print) > print_every_sec || k == Nraw
            print_progress(k, Nraw, t_now);
            t_last_print = t_now;
        end
        continue;
    end

    key = make_param_key4(m.T, m.doping, m.mu, pol);
    this_time = files(k).datenum;

    if ~isKey(meta_map, key)
        meta_map(key) = struct( ...
            "path",   p, ...
            "T",      m.T, ...
            "doping", m.doping, ...
            "mu",     m.mu, ...
            "polar",  pol, ...
            "mode",   m.mode, ...
            "time",   this_time);
    else
        old = meta_map(key);
        if this_time < old.time
            meta_map(key) = struct( ...
                "path",   p, ...
                "T",      m.T, ...
                "doping", m.doping, ...
                "mu",     m.mu, ...
                "polar",  pol, ...
                "mode",   m.mode, ...
                "time",   this_time);
        end
    end

    t_now = toc(t0);
    if (t_now - t_last_print) > print_every_sec || k == Nraw
        print_progress(k, Nraw, t_now);
        t_last_print = t_now;
    end
end

fprintf("[scan done] kept=%d | skipped(no T)=%d | time=%.1fs\n", ...
    meta_map.Count, n_skip_noT, toc(t0));

% -----------------------------
% map -> struct array meta
% -----------------------------
keys_list = meta_map.keys;
meta = repmat(struct("path","", "T",NaN, "doping",NaN, "mu",NaN, "polar",NaN, "mode",""), numel(keys_list), 1);

for i = 1:numel(keys_list)
    tmp = meta_map(keys_list{i});
    meta(i).path   = tmp.path;
    meta(i).T      = tmp.T;
    meta(i).doping = tmp.doping;
    meta(i).mu     = tmp.mu;
    meta(i).polar  = tmp.polar;
    meta(i).mode   = tmp.mode;
end

% -----------------------------
% unique T / polar
% -----------------------------
Tall = [meta.T]';
Tuniq = unique(Tall(isfinite(Tall)));
Tuniq = sort(Tuniq);
if isempty(Tuniq), error("No finite T parsed from headers."); end

Pall = [meta.polar]';
Puniq = unique(Pall(isfinite(Pall)));
Puniq = sort(Puniq);
if isempty(Puniq)
    warning("No polar parsed from path. Use polar=0 as fallback.");
    Puniq = 0;
    for i = 1:numel(meta)
        meta(i).polar = 0;
    end
end

it0 = 1;
ip0 = 1;

% -----------------------------
% plot defaults
% -----------------------------
FS = 16;
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
set(groot, "defaultTextInterpreter", "latex");

% -----------------------------
% build figure
% -----------------------------
fig = figure("Color","w", "Units","pixels", "Position",[80 80 1100 720], ...
    "Name","chi vs mu (T + polar sliders)");
ax = axes(fig, "Position",[0.10 0.26 0.86 0.66]);

hInfo = uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.945 0.86 0.035], "String","", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");

% -----------------------------
% controls
% -----------------------------
uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.16 0.08 0.035], "String","T idx", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
sT = uicontrol(fig, "Style","slider", "Units","normalized", ...
    "Position",[0.18 0.165 0.38 0.030], ...
    "Min",1, "Max",numel(Tuniq), "Value", it0, ...
    "SliderStep", slider_step(numel(Tuniq)), ...
    "Callback", @onTChanged);
tT = uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.58 0.160 0.12 0.035], "String","", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");

uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.11 0.08 0.035], "String","polar idx", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
sP = uicontrol(fig, "Style","slider", "Units","normalized", ...
    "Position",[0.18 0.115 0.38 0.030], ...
    "Min",1, "Max",numel(Puniq), "Value", ip0, ...
    "SliderStep", slider_step(numel(Puniq)), ...
    "Callback", @onPolarChanged);
tP = uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.58 0.110 0.12 0.035], "String","", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");

uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.75 0.16 0.05 0.035], "String","iq", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
eIq = uicontrol(fig, "Style","edit", "Units","normalized", ...
    "Position",[0.79 0.165 0.06 0.032], "String", num2str(iq0), ...
    "FontSize", 11, "Callback", @onIqJqChanged);

uicontrol(fig, "Style","text", "Units","normalized", ...
    "Position",[0.86 0.16 0.05 0.035], "String","jq", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
eJq = uicontrol(fig, "Style","edit", "Units","normalized", ...
    "Position",[0.90 0.165 0.06 0.032], "String", num2str(jq0), ...
    "FontSize", 11, "Callback", @onIqJqChanged);

uicontrol(fig, "Style","pushbutton", "Units","normalized", ...
    "Position",[0.78 0.060 0.18 0.050], ...
    "String","Save PNG", "FontSize", 12, ...
    "Callback", @onSave);

% -----------------------------
% shared state
% -----------------------------
state = struct();
state.meta  = meta;
state.Tuniq = Tuniq;
state.Puniq = Puniq;
state.it    = it0;
state.ip    = ip0;
state.iq    = iq0;
state.jq    = jq0;

% initial draw
update_figure();

% ============================================================
% callbacks
% ============================================================
function onTChanged(~,~)
    state.it = clampi(round(get(sT,"Value")), 1, numel(state.Tuniq));
    set(sT,"Value", state.it);
    update_figure();
end

function onPolarChanged(~,~)
    state.ip = clampi(round(get(sP,"Value")), 1, numel(state.Puniq));
    set(sP,"Value", state.ip);
    update_figure();
end

function onIqJqChanged(~,~)
    [iq, jq] = read_iqjq(eIq, eJq);
    state.iq = iq;
    state.jq = jq;
    update_figure();
end

function onSave(~,~)
    Tnow = state.Tuniq(state.it);
    Pnow = state.Puniq(state.ip);
    out_png = fullfile(root, sprintf("chi_vs_mu_T%.6g_polar%.6g_iq%d_jq%d.png", ...
        Tnow, Pnow, state.iq, state.jq));
    exportgraphics(fig, out_png, "Resolution", 300);
    fprintf("Saved:\n  %s\n", out_png);
end

% ============================================================
% core update
% ============================================================
function update_figure()
    Tnow = state.Tuniq(state.it);
    Pnow = state.Puniq(state.ip);

    set(tT,"String", sprintf("T=%.6g", Tnow));
    set(tP,"String", sprintf("polar=%.6g", Pnow));

    Ttol = 1e-9;
    Ptol = 1e-12;

    idxTP = find(abs([state.meta.T] - Tnow) < Ttol & abs([state.meta.polar] - Pnow) < Ptol);

    if isempty(idxTP)
        cla(ax);
        text(ax, 0.5, 0.5, "No files at this (T, polar)", ...
            "Units","normalized", "HorizontalAlignment","center");
        set(hInfo, "String", sprintf("T=%.12g | polar=%.12g | iq=%d jq=%d | valid=0", ...
            Tnow, Pnow, state.iq, state.jq));
        return;
    end

    paths = string({state.meta(idxTP).path});
    mu    = [state.meta(idxTP).mu]';
    chi   = nan(numel(paths),1);

    for k2 = 1:numel(paths)
        chi(k2) = extract_chi_at_iqjq(paths(k2), state.iq, state.jq);
    end

    ok = isfinite(chi) & isfinite(mu);
    paths = paths(ok); %#ok<NASGU>
    mu    = mu(ok);
    chi   = chi(ok);

    set(hInfo, "String", sprintf("T=%.12g | polar=%.12g | iq=%d jq=%d | valid=%d | plot: \\chi(\\mu)", ...
        Tnow, Pnow, state.iq, state.jq, numel(chi)));

    cla(ax);

    if isempty(chi)
        text(ax, 0.5, 0.5, "No valid chi or mu at this (T, polar, iq, jq)", ...
            "Units","normalized", "HorizontalAlignment","center");
        return;
    end

    hold(ax,"on");

    scatter(ax, mu, chi, 18, "filled", ...
        "MarkerFaceAlpha", 0.20, "MarkerEdgeAlpha", 0.20);

    [xm, ym_mean, ~, ~] = group_mean_std(mu, chi);
    if ~isempty(xm)
        plot(ax, xm, ym_mean, "-", "LineWidth", 1.8);
    end

    hold(ax,"off");

    xlabel(ax, "\mu (eV)", "FontSize", FS);
    ylabel(ax, "Re \chi(iq,jq)", "FontSize", FS);
    title(ax, sprintf("Re \\chi vs \\mu @ T=%.6g, polar=%.6g (iq=%d,jq=%d)", ...
        Tnow, Pnow, state.iq, state.jq), ...
        "FontSize", FS, "FontWeight","normal");
    grid(ax,"on");
    set(ax, "FontSize", FS, "LineWidth", 1.0, "TickDir","out", "Box","on");
end

end % ===== end main =====


% ============================================================
% helpers
% ============================================================

function pol = parse_polar_from_path(p)
pol = NaN;
s = char(p);
tok = regexp(s, "polar[_\-]meV([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)", "tokens", "once");
if ~isempty(tok)
    pol = str2double(tok{1});
end
end

function [chi_run_dir, Tlist] = resolve_T_root(D_dir)
Tlist = dir(fullfile(D_dir, "T*"));
Tlist = Tlist([Tlist.isdir]);
Tlist = Tlist(~ismember({Tlist.name},{'.','..'}));
if ~isempty(Tlist)
    chi_run_dir = string(D_dir);
    return;
end

sub = dir(D_dir);
sub = sub([sub.isdir]);
sub = sub(~ismember({sub.name},{'.','..'}));

for i = 1:numel(sub)
    cand = fullfile(D_dir, sub(i).name);
    T2 = dir(fullfile(cand, "T*"));
    T2 = T2([T2.isdir]);
    T2 = T2(~ismember({T2.name},{'.','..'}));
    if ~isempty(T2)
        chi_run_dir = string(D_dir);
        Tlist = T2;
        return;
    end
end

error("Cannot find any T* folders under:\n  %s\nor under one subfolder level.", D_dir);
end

function step = slider_step(n)
if n <= 1
    step = [1 1];
else
    step = [1/(n-1) min(10/(n-1),1)];
end
end

function x = clampi(x, a, b)
x = min(max(x, a), b);
end

function [iq, jq] = read_iqjq(eIq, eJq)
iq = str2double(get(eIq,"String"));
jq = str2double(get(eJq,"String"));
if ~isfinite(iq), iq = 0; end
if ~isfinite(jq), jq = 0; end
iq = round(iq);
jq = round(jq);
set(eIq,"String", num2str(iq));
set(eJq,"String", num2str(jq));
end

function key = make_param_key4(T, dop, mu, pol)
key = sprintf("T=%s|dop=%s|mu=%s|pol=%s", fnum(T), fnum(dop), fnum(mu), fnum(pol));
end

function s = fnum(x)
if ~isfinite(x)
    s = "nan";
else
    s = sprintf("%.12g", x);
end
end

function print_progress(k, N, elapsed)
pct = 100*k/N;
if k <= 1
    eta = NaN;
else
    rate = elapsed / k;
    eta  = rate * (N - k);
end
if isfinite(eta)
    fprintf("  scan %4d/%4d (%.1f%%)  elapsed %.1fs  ETA %.1fs\n", k, N, pct, elapsed, eta);
else
    fprintf("  scan %4d/%4d (%.1f%%)  elapsed %.1fs\n", k, N, pct, elapsed);
end
drawnow;
end

function m = parse_header_T_doping_mu_mode(in_path)
m = struct("T",NaN,"doping",NaN,"mu",NaN,"mode","");

fid = fopen(in_path,'r');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

while true
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    s0 = strtrim(tline);
    if ~startsWith(s0, "#")
        break;
    end
    s = strtrim(erase(s0, "#"));

    if strlength(m.mode)==0
        tok = regexp(s, "mode\s*=\s*([A-Za-z0-9_\-]+)", "tokens", "once");
        if ~isempty(tok), m.mode = string(tok{1}); end
    end

    if ~isfinite(m.T)
        tok = regexp(s, "(?:^|\s)T(?:_K)?\s*=\s*" + num, "tokens", "once");
        if ~isempty(tok), m.T = str2double(tok{1}); end
    end

    if ~isfinite(m.doping)
        tok = regexp(s, "doping\s*=\s*" + num, "tokens", "once");
        if ~isempty(tok), m.doping = str2double(tok{1}); end
    end

    if ~isfinite(m.mu)
        tok = regexp(s, "(?:^|\s)(?:mu|mu_eV|EF|EF_eV)\s*=\s*" + num, "tokens", "once");
        if ~isempty(tok), m.mu = str2double(tok{1}); end
    end
end
end

function val = extract_chi_at_iqjq(in_path, iq0, jq0)
val = NaN;

try
    raw = readmatrix(in_path, "FileType","text", "CommentStyle","#");
catch
    return;
end

if isempty(raw) || size(raw,2) < 6
    return;
end

if size(raw,2) >= 8
    iq = raw(:,2);
    jq = raw(:,3);
    re = raw(:,6);
else
    iq = raw(:,1);
    jq = raw(:,2);
    re = raw(:,5);
end

idx = find(iq == iq0 & jq == jq0, 1, "first");
if isempty(idx), return; end
val = double(re(idx));
end

function [xuniq, ymean, ystd, ncount] = group_mean_std(x, y)
x = double(x(:));
y = double(y(:));

ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);

if isempty(x)
    xuniq = [];
    ymean = [];
    ystd = [];
    ncount = [];
    return;
end

[xuniq, ~, ic] = unique(x);
ymean  = accumarray(ic, y, [], @mean);
ystd   = accumarray(ic, y, [], @std);
ncount = accumarray(ic, 1, [], @sum);

[xuniq, ord] = sort(xuniq);
ymean  = ymean(ord);
ystd   = ystd(ord); 
ncount = ncount(ord);
end