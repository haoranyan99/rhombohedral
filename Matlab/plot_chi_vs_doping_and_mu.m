function out_png = plot_chi_vs_doping_and_mu()
% plot_chi_vs_doping_and_mu2
% Scan chi*.txt recursively under a D folder, parse header T/doping/mu(EF),
% deduplicate by (T,doping,mu) keeping the OLDEST file, then build TWO figures:
%   - Figure D: Re chi(iq,jq) vs doping  (T slider)
%   - Figure M: Re chi(iq,jq) vs mu      (T slider)
%
% Robust to mixed run dates and inconsistent filenames: we only rely on headers
% and file datenum (oldest kept for each parameter set).
%
% File numeric table formats supported:
%   if >=8 cols: idx iq jq qx qy Re Im ...
%   else:        iq jq qx qy Re Im ...
%
% NEW: user-defined filter by mu OR doping range; skip files outside range.

out_png = "";

% =============================
% USER FILTER (EDIT HERE)
% =============================
filter = struct();
filter.mode  = "mu";          % "mu" or "doping"
filter.range = [1.0045, 2]; % [min max] (inclusive)
% =============================

% -----------------------------
% default iq/jq
% -----------------------------
iq0 = -133;
jq0 = -133;

% -----------------------------
% choose root folder
% -----------------------------
default_root = "/Users/haoranyan/Library/CloudStorage/OneDrive-Emory/Rhombohedral_SC/rhombohedral_project/data/chi_sk_mu_200/D0.067";
default_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_200\D0.067";
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
fprintf("[filter] mode=%s | range=[%.12g, %.12g]\n", filter.mode, filter.range(1), filter.range(2));

% -----------------------------
% parse headers and keep OLDEST per (T,doping,mu)
% -----------------------------
meta_map = containers.Map('KeyType','char','ValueType','any');

t0 = tic;
t_last_print = 0;
print_every_sec = 0.5; % update progress about twice per sec

n_used = 0;
n_skip_noT = 0;
n_skip_filter = 0;

for k = 1:Nraw
    p = string(fullfile(files(k).folder, files(k).name));

    m = parse_header_T_doping_mu(p);

    if ~isfinite(m.T)
        n_skip_noT = n_skip_noT + 1;
        t_now = toc(t0);
        if (t_now - t_last_print) > print_every_sec || k==Nraw
            print_progress(k, Nraw, t_now);
            t_last_print = t_now;
        end
        continue;
    end

    % ---------- NEW: range filter (early skip) ----------
    if ~pass_filter(m, filter)
        n_skip_filter = n_skip_filter + 1;
        t_now = toc(t0);
        if (t_now - t_last_print) > print_every_sec || k==Nraw
            print_progress(k, Nraw, t_now);
            t_last_print = t_now;
        end
        continue;
    end
    % ---------------------------------------------------

    % key: T, doping, mu (can be NaN; keep it as NaN string)
    key = make_param_key(m.T, m.doping, m.mu);

    this_time = files(k).datenum;

    if ~isKey(meta_map, key)
        meta_map(key) = struct("path", p, "T", m.T, "doping", m.doping, "mu", m.mu, "time", this_time);
        n_used = n_used + 1;
    else
        old = meta_map(key);
        if this_time < old.time
            meta_map(key) = struct("path", p, "T", m.T, "doping", m.doping, "mu", m.mu, "time", this_time);
        end
    end

    % progress printing
    t_now = toc(t0);
    if (t_now - t_last_print) > print_every_sec || k==Nraw
        print_progress(k, Nraw, t_now);
        t_last_print = t_now;
    end
end

fprintf("[scan done] kept=%d | skipped(no T)=%d | skipped(filter)=%d | time=%.1fs\n", ...
    meta_map.Count, n_skip_noT, n_skip_filter, toc(t0));

% -----------------------------
% map -> struct array meta
% -----------------------------
keys_list = meta_map.keys;
meta = repmat(struct("path","", "T",NaN, "doping",NaN, "mu",NaN), numel(keys_list), 1);

for i = 1:numel(keys_list)
    tmp = meta_map(keys_list{i});
    meta(i).path   = tmp.path;
    meta(i).T      = tmp.T;
    meta(i).doping = tmp.doping;
    meta(i).mu     = tmp.mu;
end

% -----------------------------
% unique T from parsed headers
% -----------------------------
Tall = [meta.T]';
Tuniq = unique(Tall(isfinite(Tall)));
Tuniq = sort(Tuniq);
if isempty(Tuniq)
    error("No finite T parsed from headers (maybe all filtered out?).");
end
it0 = 1;

% -----------------------------
% plot defaults
% -----------------------------
FS = 16;
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
set(groot, "defaultTextInterpreter", "latex");

% -----------------------------
% build TWO figures
% -----------------------------
figD = figure("Color","w", "Units","pixels", "Position",[80 80 1100 720], ...
    "Name","chi vs doping (T slider)");
axD = axes(figD, "Position",[0.10 0.22 0.86 0.70]);

figM = figure("Color","w", "Units","pixels", "Position",[120 120 1100 720], ...
    "Name","chi vs mu (T slider)");
axM = axes(figM, "Position",[0.10 0.22 0.86 0.70]);

% top info text
hInfoD = uicontrol(figD, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.945 0.86 0.035], "String","", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
hInfoM = uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.945 0.86 0.035], "String","", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");

% -----------------------------
% controls (Figure D)
% -----------------------------
uicontrol(figD, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.13 0.08 0.035], "String","T idx", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
sTD = uicontrol(figD, "Style","slider", "Units","normalized", ...
    "Position",[0.18 0.135 0.62 0.030], ...
    "Min",1, "Max",numel(Tuniq), "Value", it0, ...
    "SliderStep", slider_step(numel(Tuniq)), ...
    "Callback", @onTChanged_from_D);

uicontrol(figD, "Style","text", "Units","normalized", ...
    "Position",[0.82 0.13 0.05 0.035], "String","iq", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
eIqD = uicontrol(figD, "Style","edit", "Units","normalized", ...
    "Position",[0.86 0.135 0.05 0.032], "String", num2str(iq0), ...
    "FontSize", 11, "Callback", @onIqJqChanged_from_D);

uicontrol(figD, "Style","text", "Units","normalized", ...
    "Position",[0.92 0.13 0.05 0.035], "String","jq", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
eJqD = uicontrol(figD, "Style","edit", "Units","normalized", ...
    "Position",[0.96 0.135 0.05 0.032], "String", num2str(jq0), ...
    "FontSize", 11, "Callback", @onIqJqChanged_from_D);

% -----------------------------
% controls (Figure M)
% -----------------------------
uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.13 0.08 0.035], "String","T idx", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
sTM = uicontrol(figM, "Style","slider", "Units","normalized", ...
    "Position",[0.18 0.135 0.62 0.030], ...
    "Min",1, "Max",numel(Tuniq), "Value", it0, ...
    "SliderStep", slider_step(numel(Tuniq)), ...
    "Callback", @onTChanged_from_M);

uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.82 0.13 0.05 0.035], "String","iq", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
eIqM = uicontrol(figM, "Style","edit", "Units","normalized", ...
    "Position",[0.86 0.135 0.05 0.032], "String", num2str(iq0), ...
    "FontSize", 11, "Callback", @onIqJqChanged_from_M);

uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.92 0.13 0.05 0.035], "String","jq", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
eJqM = uicontrol(figM, "Style","edit", "Units","normalized", ...
    "Position",[0.96 0.135 0.05 0.032], "String", num2str(jq0), ...
    "FontSize", 11, "Callback", @onIqJqChanged_from_M);

% -----------------------------
% shared state
% -----------------------------
state = struct();
state.meta   = meta;
state.Tuniq  = Tuniq;
state.it     = it0;
state.iq     = iq0;
state.jq     = jq0;
state.filter = filter;

% initial draw
update_both_figures();

% ============================================================
% callbacks
% ============================================================
function onTChanged_from_D(~,~)
    state.it = clampi(round(get(sTD,"Value")), 1, numel(state.Tuniq));
    set(sTD,"Value", state.it);
    set(sTM,"Value", state.it);
    update_both_figures();
end

function onTChanged_from_M(~,~)
    state.it = clampi(round(get(sTM,"Value")), 1, numel(state.Tuniq));
    set(sTM,"Value", state.it);
    set(sTD,"Value", state.it);
    update_both_figures();
end

function onIqJqChanged_from_D(~,~)
    [iq, jq] = read_iqjq(eIqD, eJqD);
    state.iq = iq; state.jq = jq;
    set(eIqM,"String", num2str(iq));
    set(eJqM,"String", num2str(jq));
    update_both_figures();
end

function onIqJqChanged_from_M(~,~)
    [iq, jq] = read_iqjq(eIqM, eJqM);
    state.iq = iq; state.jq = jq;
    set(eIqD,"String", num2str(iq));
    set(eJqD,"String", num2str(jq));
    update_both_figures();
end

% ============================================================
% core update
% ============================================================
function update_both_figures()
    Tnow = state.Tuniq(state.it);

    % exact match (header generated, should match exactly; allow tiny tol)
    Ttol = 1e-9;
    idxT = find(abs([state.meta.T] - Tnow) < Ttol);

    if isempty(idxT)
        % fallback: nearest T
        [~, j] = min(abs([state.meta.T] - Tnow));
        Tnow = state.meta(j).T;
        idxT = find(abs([state.meta.T] - Tnow) < Ttol);
    end

    paths = string({state.meta(idxT).path});
    dop   = [state.meta(idxT).doping]';
    mu    = [state.meta(idxT).mu]';

    % extract chi at (iq,jq)
    chi = nan(numel(paths),1);
    for k2 = 1:numel(paths)
        chi(k2) = extract_chi_at_iqjq(paths(k2), state.iq, state.jq);
    end

    ok = isfinite(chi);
    paths = paths(ok);
    dop   = dop(ok);
    mu    = mu(ok);
    chi   = chi(ok);

    fstr = sprintf("filter: %s in [%.6g, %.6g]", state.filter.mode, state.filter.range(1), state.filter.range(2));
    set(hInfoD, "String", sprintf("T=%.12g K | iq=%d jq=%d | valid files=%d | %s", Tnow, state.iq, state.jq, numel(chi), fstr));
    set(hInfoM, "String", sprintf("T=%.12g K | iq=%d jq=%d | valid files=%d | %s", Tnow, state.iq, state.jq, numel(chi), fstr));

    % -------- Figure 1: chi vs doping --------
    cla(axD);
    if isempty(chi)
        text(axD, 0.5, 0.5, "No valid chi at this (T,iq,jq) under filter", "Units","normalized", ...
            "HorizontalAlignment","center");
    else
        hold(axD,"on");
        if any(isfinite(dop))
            scatter(axD, dop, chi, 18, "filled", "MarkerFaceAlpha", 0.25, "MarkerEdgeAlpha", 0.25);
        end

        okd = isfinite(dop);
        [xd, yd_mean, ~, ~] = group_mean_std(dop(okd), chi(okd));
        if ~isempty(xd)
            plot(axD, xd, yd_mean, "-", "LineWidth", 1.8);
        end
        hold(axD,"off");

        xlabel(axD, "doping (10^{12} cm^{-2})", "FontSize", FS);
        ylabel(axD, "Re \chi(iq,jq)", "FontSize", FS);
        title(axD, sprintf("Re \\chi vs doping @ T=%.6g K (iq=%d, jq=%d)", Tnow, state.iq, state.jq), ...
            "FontSize", FS, "FontWeight","normal");
        grid(axD,"on");
        set(axD, "FontSize", FS, "LineWidth", 1.0, "TickDir","out", "Box","on");
    end

    % -------- Figure 2: chi vs mu --------
    cla(axM);
    if isempty(chi)
        text(axM, 0.5, 0.5, "No valid chi at this (T,iq,jq) under filter", "Units","normalized", ...
            "HorizontalAlignment","center");
    else
        hold(axM,"on");
        if any(isfinite(mu))
            scatter(axM, mu, chi, 18, "filled", "MarkerFaceAlpha", 0.25, "MarkerEdgeAlpha", 0.25);
        end

        okm = isfinite(mu);
        [xm, ym_mean, ~, ~] = group_mean_std(mu(okm), chi(okm));
        if ~isempty(xm)
            plot(axM, xm, ym_mean, "-", "LineWidth", 1.8);
        end
        hold(axM,"off");

        xlabel(axM, "\mu (eV)", "FontSize", FS);
        ylabel(axM, "Re \chi(iq,jq)", "FontSize", FS);
        title(axM, sprintf("Re \\chi vs \\mu @ T=%.6g K (iq=%d, jq=%d)", Tnow, state.iq, state.jq), ...
            "FontSize", FS, "FontWeight","normal");
        grid(axM,"on");
        set(axM, "FontSize", FS, "LineWidth", 1.0, "TickDir","out", "Box","on");
    end
end

end % ===== end main =====


% ============================================================
% helpers (ALL included)
% ============================================================

function tf = pass_filter(m, filter)
    % Return true if this file passes the (mu OR doping) range filter.
    lo = min(filter.range); hi = max(filter.range);

    mode = lower(string(filter.mode));
    if mode == "mu"
        v = m.mu;
    elseif mode == "doping"
        v = m.doping;
    else
        error('filter.mode must be "mu" or "doping".');
    end

    if ~isfinite(v)
        tf = false; % missing value -> skip (range filtering requires it)
        return;
    end
    tf = (v >= lo) && (v <= hi);
end

function [chi_run_dir, Tlist] = resolve_T_root(D_dir)
    % Case A: directly contains T*
    Tlist = dir(fullfile(D_dir, "T*"));
    Tlist = Tlist([Tlist.isdir]);
    Tlist = Tlist(~ismember({Tlist.name},{'.','..'}));
    if ~isempty(Tlist)
        chi_run_dir = string(D_dir);
        return;
    end

    % Case B: one extra run_tag layer
    sub = dir(D_dir);
    sub = sub([sub.isdir]);
    sub = sub(~ismember({sub.name},{'.','..'}));

    for i = 1:numel(sub)
        cand = fullfile(D_dir, sub(i).name);
        T2 = dir(fullfile(cand, "T*"));
        T2 = T2([T2.isdir]);
        T2 = T2(~ismember({T2.name},{'.','..'}));
        if ~isempty(T2)
            chi_run_dir = string(cand);
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

function key = make_param_key(T, dop, mu)
    % stable string key; keeps NaN as "nan"
    key = sprintf("T=%s|dop=%s|mu=%s", fnum(T), fnum(dop), fnum(mu));
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
        rate = elapsed / k;           % sec per file
        eta  = rate * (N - k);
    end
    if isfinite(eta)
        fprintf("  scan %4d/%4d (%.1f%%)  elapsed %.1fs  ETA %.1fs\n", k, N, pct, elapsed, eta);
    else
        fprintf("  scan %4d/%4d (%.1f%%)  elapsed %.1fs\n", k, N, pct, elapsed);
    end
    drawnow;
end

function m = parse_header_T_doping_mu(in_path)
    % Parse header lines starting with '#'
    % Accept:
    %   T =, T_K =
    %   doping =
    %   mu =, mu_eV =, EF =, EF_eV =
    m = struct("T",NaN,"doping",NaN,"mu",NaN);

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
    % Return Re(chi) at target (iq,jq), or NaN if missing.
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
    x = x(ok); y = y(ok);
    if isempty(x)
        xuniq = []; ymean = []; ystd = []; ncount = [];
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
