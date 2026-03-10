function out_png = plot_chi_vs_doping_and_mu()
% plot_chi_vs_doping_and_mu (with T + polar + mu sliders)
% Scan chi*.txt recursively under a D folder, parse header T/doping/mu(EF),
% also parse "polar" from folder name like: polar_meV0.000
% ALSO parse "mode" from header line like: "# mode = mu"
% Deduplicate by (T,doping,mu,polar) keeping the OLDEST file.
%
% Build TWO figures:
%   - Figure D: Re chi(iq,jq) vs doping  (T slider + polar slider + mu slider)
%       mu slider FIXES one mu (nearest available) and shows chi(doping) at that mu.
%   - Figure M:
%       if mode != mu: Re chi(iq,jq) vs mu      (T slider + polar slider + mu slider marker)
%       if mode == mu: Re chi(iq,jq) vs doping  (T slider + polar slider)  [use header doping]
%
% File numeric table formats supported:
%   if >=8 cols: idx iq jq qx qy Re Im ...
%   else:        iq jq qx qy Re Im ...

out_png = "";

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

% -----------------------------
% parse headers and keep OLDEST per (T,doping,mu,polar)
% also collect mode strings
% -----------------------------
meta_map = containers.Map('KeyType','char','ValueType','any');
mode_list = strings(0,1);

t0 = tic;
t_last_print = 0;
print_every_sec = 0.5;

n_used = 0;
n_skip_noT = 0;

for k = 1:Nraw
    p = string(fullfile(files(k).folder, files(k).name));

    m = parse_header_T_doping_mu_mode(p);
    pol = parse_polar_from_path(p);

    if ~isfinite(m.T)
        n_skip_noT = n_skip_noT + 1;
        t_now = toc(t0);
        if (t_now - t_last_print) > print_every_sec || k==Nraw
            print_progress(k, Nraw, t_now);
            t_last_print = t_now;
        end
        continue;
    end

    if strlength(m.mode) > 0
        mode_list(end+1,1) = lower(strtrim(m.mode)); %#ok<AGROW>
    end

    key = make_param_key4(m.T, m.doping, m.mu, pol);
    this_time = files(k).datenum;

    if ~isKey(meta_map, key)
        meta_map(key) = struct("path", p, "T", m.T, "doping", m.doping, ...
            "mu", m.mu, "polar", pol, "mode", m.mode, "time", this_time);
        n_used = n_used + 1;
    else
        old = meta_map(key);
        if this_time < old.time
            meta_map(key) = struct("path", p, "T", m.T, "doping", m.doping, ...
                "mu", m.mu, "polar", pol, "mode", m.mode, "time", this_time);
        end
    end

    t_now = toc(t0);
    if (t_now - t_last_print) > print_every_sec || k==Nraw
        print_progress(k, Nraw, t_now);
        t_last_print = t_now;
    end
end

fprintf("[scan done] kept=%d | skipped(no T)=%d | time=%.1fs\n", ...
    meta_map.Count, n_skip_noT, toc(t0));

% pick global mode (most frequent); fallback ""
mode_global = "";
if ~isempty(mode_list)
    [u, ~, ic] = unique(mode_list);
    cnt = accumarray(ic, 1);
    [~, imax] = max(cnt);
    mode_global = u(imax);
end
fprintf("[mode] global mode detected = '%s'\n", mode_global);

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
% unique T / polar / mu
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
    for i=1:numel(meta), meta(i).polar = 0; end
end

Muall = [meta.mu]';
Muuniq = unique(Muall(isfinite(Muall)));
Muuniq = sort(Muuniq);
if isempty(Muuniq)
    warning("No mu/EF parsed from headers; mu slider will be disabled.");
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
% build TWO figures
% -----------------------------
figD = figure("Color","w", "Units","pixels", "Position",[80 80 1100 720], ...
    "Name","chi vs doping (T+polar+mu sliders)");
axD = axes(figD, "Position",[0.10 0.26 0.86 0.66]);

figM = figure("Color","w", "Units","pixels", "Position",[120 120 1100 720], ...
    "Name","chi (mode-adaptive) (T+polar+mu sliders)");
axM = axes(figM, "Position",[0.10 0.26 0.86 0.66]);

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
    "Position",[0.10 0.16 0.08 0.035], "String","T idx", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
sTD = uicontrol(figD, "Style","slider", "Units","normalized", ...
    "Position",[0.18 0.165 0.38 0.030], ...
    "Min",1, "Max",numel(Tuniq), "Value", it0, ...
    "SliderStep", slider_step(numel(Tuniq)), ...
    "Callback", @onTChanged_from_D);
tTD = uicontrol(figD, "Style","text", "Units","normalized", ...
    "Position",[0.58 0.160 0.12 0.035], "String","", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");

uicontrol(figD, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.11 0.08 0.035], "String","polar idx", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
sPD = uicontrol(figD, "Style","slider", "Units","normalized", ...
    "Position",[0.18 0.115 0.38 0.030], ...
    "Min",1, "Max",numel(Puniq), "Value", ip0, ...
    "SliderStep", slider_step(numel(Puniq)), ...
    "Callback", @onPolarChanged_from_D);
tPD = uicontrol(figD, "Style","text", "Units","normalized", ...
    "Position",[0.58 0.110 0.12 0.035], "String","", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");

uicontrol(figD, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.06 0.08 0.035], "String","mu idx", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
sMuD = uicontrol(figD, "Style","slider", "Units","normalized", ...
    "Position",[0.18 0.065 0.38 0.030], ...
    "Min",1, "Max",max(numel(Muuniq),1), "Value", 1, ...
    "SliderStep", slider_step(max(numel(Muuniq),1)), ...
    "Callback", @onMuChanged_from_D);
tMuD = uicontrol(figD, "Style","text", "Units","normalized", ...
    "Position",[0.58 0.060 0.12 0.035], "String","", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");

uicontrol(figD, "Style","text", "Units","normalized", ...
    "Position",[0.75 0.16 0.05 0.035], "String","iq", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
eIqD = uicontrol(figD, "Style","edit", "Units","normalized", ...
    "Position",[0.79 0.165 0.06 0.032], "String", num2str(iq0), ...
    "FontSize", 11, "Callback", @onIqJqChanged_from_D);

uicontrol(figD, "Style","text", "Units","normalized", ...
    "Position",[0.86 0.16 0.05 0.035], "String","jq", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
eJqD = uicontrol(figD, "Style","edit", "Units","normalized", ...
    "Position",[0.90 0.165 0.06 0.032], "String", num2str(jq0), ...
    "FontSize", 11, "Callback", @onIqJqChanged_from_D);

% -----------------------------
% controls (Figure M)
% -----------------------------
uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.16 0.08 0.035], "String","T idx", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
sTM = uicontrol(figM, "Style","slider", "Units","normalized", ...
    "Position",[0.18 0.165 0.38 0.030], ...
    "Min",1, "Max",numel(Tuniq), "Value", it0, ...
    "SliderStep", slider_step(numel(Tuniq)), ...
    "Callback", @onTChanged_from_M);
tTM = uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.58 0.160 0.12 0.035], "String","", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");

uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.11 0.08 0.035], "String","polar idx", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
sPM = uicontrol(figM, "Style","slider", "Units","normalized", ...
    "Position",[0.18 0.115 0.38 0.030], ...
    "Min",1, "Max",numel(Puniq), "Value", ip0, ...
    "SliderStep", slider_step(numel(Puniq)), ...
    "Callback", @onPolarChanged_from_M);
tPM = uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.58 0.110 0.12 0.035], "String","", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");

uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.10 0.06 0.08 0.035], "String","mu idx", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
sMuM = uicontrol(figM, "Style","slider", "Units","normalized", ...
    "Position",[0.18 0.065 0.38 0.030], ...
    "Min",1, "Max",max(numel(Muuniq),1), "Value", 1, ...
    "SliderStep", slider_step(max(numel(Muuniq),1)), ...
    "Callback", @onMuChanged_from_M);
tMuM = uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.58 0.060 0.12 0.035], "String","", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");

uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.75 0.16 0.05 0.035], "String","iq", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
eIqM = uicontrol(figM, "Style","edit", "Units","normalized", ...
    "Position",[0.79 0.165 0.06 0.032], "String", num2str(iq0), ...
    "FontSize", 11, "Callback", @onIqJqChanged_from_M);

uicontrol(figM, "Style","text", "Units","normalized", ...
    "Position",[0.86 0.16 0.05 0.035], "String","jq", ...
    "BackgroundColor","w", "FontSize", 11, "HorizontalAlignment","left");
eJqM = uicontrol(figM, "Style","edit", "Units","normalized", ...
    "Position",[0.90 0.165 0.06 0.032], "String", num2str(jq0), ...
    "FontSize", 11, "Callback", @onIqJqChanged_from_M);

% disable mu sliders if no mu parsed
if isempty(Muuniq)
    set(sMuD,"Enable","off"); set(sMuM,"Enable","off");
    set(tMuD,"String","(no mu)"); set(tMuM,"String","(no mu)");
end

% if mode==mu, Figure M is vs doping -> mu slider is not essential; we disable it to avoid confusion
if strcmpi(mode_global, "mu")
    set(sMuM, "Enable", "off");
    set(tMuM, "String", "(mode=mu: figM plots vs doping)");
end

% -----------------------------
% shared state
% -----------------------------
state = struct();
state.meta   = meta;
state.Tuniq  = Tuniq;
state.Puniq  = Puniq;
state.Muuniq = Muuniq;
state.mode_global = mode_global;

state.it     = it0;
state.ip     = ip0;
state.imu    = 1;

state.iq     = iq0;
state.jq     = jq0;

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

function onPolarChanged_from_D(~,~)
    state.ip = clampi(round(get(sPD,"Value")), 1, numel(state.Puniq));
    set(sPD,"Value", state.ip);
    set(sPM,"Value", state.ip);
    update_both_figures();
end

function onPolarChanged_from_M(~,~)
    state.ip = clampi(round(get(sPM,"Value")), 1, numel(state.Puniq));
    set(sPM,"Value", state.ip);
    set(sPD,"Value", state.ip);
    update_both_figures();
end

function onMuChanged_from_D(~,~)
    if isempty(state.Muuniq), return; end
    state.imu = clampi(round(get(sMuD,"Value")), 1, numel(state.Muuniq));
    set(sMuD,"Value", state.imu);
    set(sMuM,"Value", state.imu);
    update_both_figures();
end

function onMuChanged_from_M(~,~)
    if isempty(state.Muuniq), return; end
    if strcmpi(state.mode_global, "mu")
        % figM not using mu slider in mode=mu
        return;
    end
    state.imu = clampi(round(get(sMuM,"Value")), 1, numel(state.Muuniq));
    set(sMuM,"Value", state.imu);
    set(sMuD,"Value", state.imu);
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
    Pnow = state.Puniq(state.ip);

    if isempty(state.Muuniq)
        mu_pick = NaN;
    else
        mu_pick = state.Muuniq(state.imu);
    end

    set(tTD,"String", sprintf("T=%.6g", Tnow));
    set(tTM,"String", sprintf("T=%.6g", Tnow));

    set(tPD,"String", sprintf("polar=%.6g", Pnow));
    set(tPM,"String", sprintf("polar=%.6g", Pnow));

    if ~isempty(state.Muuniq)
        set(tMuD,"String", sprintf("\\mu=%.6g", mu_pick));
        if ~strcmpi(state.mode_global, "mu")
            set(tMuM,"String", sprintf("\\mu=%.6g", mu_pick));
        end
    end

    % select meta by (T, polar)
    Ttol = 1e-9;
    Ptol = 1e-12;
    idxTP = find(abs([state.meta.T] - Tnow) < Ttol & abs([state.meta.polar] - Pnow) < Ptol);

    if isempty(idxTP)
        % fallback: nearest T, nearest polar (very simple fallback)
        [~, jT] = min(abs([state.meta.T] - Tnow));
        [~, jP] = min(abs([state.meta.polar] - Pnow));
        Tnow2 = state.meta(jT).T;
        Pnow2 = state.meta(jP).polar;
        idxTP = find(abs([state.meta.T] - Tnow2) < Ttol & abs([state.meta.polar] - Pnow2) < Ptol);
        Tnow = Tnow2; Pnow = Pnow2;
        set(tTD,"String", sprintf("T=%.6g", Tnow));
        set(tTM,"String", sprintf("T=%.6g", Tnow));
        set(tPD,"String", sprintf("polar=%.6g", Pnow));
        set(tPM,"String", sprintf("polar=%.6g", Pnow));
    end

    paths = string({state.meta(idxTP).path});
    dop   = [state.meta(idxTP).doping]';
    mu    = [state.meta(idxTP).mu]';

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

    % info
    set(hInfoD, "String", sprintf("mode='%s' | T=%.12g | polar=%.12g | iq=%d jq=%d | valid=%d | FigD: fix one \\mu then \\chi(doping)", ...
        state.mode_global, Tnow, Pnow, state.iq, state.jq, numel(chi)));

    if strcmpi(state.mode_global, "mu")
        set(hInfoM, "String", sprintf("mode='mu' | T=%.12g | polar=%.12g | iq=%d jq=%d | valid=%d | FigM: \\chi(doping) parsed from header", ...
            Tnow, Pnow, state.iq, state.jq, numel(chi)));
    else
        set(hInfoM, "String", sprintf("mode='%s' | T=%.12g | polar=%.12g | iq=%d jq=%d | valid=%d | FigM: \\chi(\\mu) and mark \\mu slider", ...
            state.mode_global, Tnow, Pnow, state.iq, state.jq, numel(chi)));
    end

    % -------- Figure D: chi vs doping at fixed mu (nearest available) --------
    cla(axD);
    if isempty(chi)
        text(axD, 0.5, 0.5, "No valid chi at this (T,polar,iq,jq)", "Units","normalized", ...
            "HorizontalAlignment","center");
    else
        if isfinite(mu_pick)
            mu_ok = mu(isfinite(mu));
            if isempty(mu_ok)
                text(axD, 0.5, 0.5, "No finite mu in valid points.", "Units","normalized", ...
                    "HorizontalAlignment","center");
            else
                % pick nearest mu among available points (BUGFIX)
                [~, jj] = min(abs(mu_ok - mu_pick));
                mu_near = mu_ok(jj);

                % use a reasonable tolerance for float comparison
                mutol = 1e-6; % was 1e-9 (too strict)
                idxMu = abs(mu - mu_near) <= mutol;

                dop2 = dop(idxMu);
                chi2 = chi(idxMu);

                if isempty(chi2)
                    text(axD, 0.5, 0.5, sprintf("No points at selected \\mu (try another mu idx). \\nUsed mutol=%.1e", mutol), ...
                        "Units","normalized", "HorizontalAlignment","center");
                else
                    hold(axD,"on");
                    scatter(axD, dop2, chi2, 22, "filled", "MarkerFaceAlpha", 0.25, "MarkerEdgeAlpha", 0.25);

                    okd = isfinite(dop2);
                    [xd, yd_mean, ~, ~] = group_mean_std(dop2(okd), chi2(okd));
                    if ~isempty(xd)
                        plot(axD, xd, yd_mean, "-", "LineWidth", 1.8);
                    end
                    hold(axD,"off");

                    xlabel(axD, "doping (10^{12} cm^{-2})", "FontSize", FS);
                    ylabel(axD, "Re \chi(iq,jq)", "FontSize", FS);
                    title(axD, sprintf("Re \\chi vs doping @ T=%.6g, polar=%.6g, \\mu\\approx%.6g (iq=%d,jq=%d)", ...
                        Tnow, Pnow, mu_near, state.iq, state.jq), "FontSize", FS, "FontWeight","normal");
                    grid(axD,"on");
                    set(axD, "FontSize", FS, "LineWidth", 1.0, "TickDir","out", "Box","on");
                end
            end
        else
            text(axD, 0.5, 0.5, "mu/EF not parsed from headers -> cannot fix mu for chi(doping)", "Units","normalized", ...
                "HorizontalAlignment","center");
        end
    end

    % -------- Figure M: mode-adaptive --------
    cla(axM);
    if isempty(chi)
        text(axM, 0.5, 0.5, "No valid chi at this (T,polar,iq,jq)", "Units","normalized", ...
            "HorizontalAlignment","center");
        return;
    end

    if strcmpi(state.mode_global, "mu")
        % === requested behavior: when mode=mu -> plot vs doping (from header) ===
        hold(axM,"on");
        okd = isfinite(dop);
        scatter(axM, dop(okd), chi(okd), 18, "filled", "MarkerFaceAlpha", 0.18, "MarkerEdgeAlpha", 0.18);
        [xd, yd_mean, ~, ~] = group_mean_std(dop(okd), chi(okd));
        if ~isempty(xd)
            plot(axM, xd, yd_mean, "-", "LineWidth", 1.8);
        end
        hold(axM,"off");

        xlabel(axM, "doping (10^{12} cm^{-2})", "FontSize", FS);
        ylabel(axM, "Re \chi(iq,jq)", "FontSize", FS);
        title(axM, sprintf("mode=mu: Re \\chi vs doping @ T=%.6g, polar=%.6g (iq=%d,jq=%d)", ...
            Tnow, Pnow, state.iq, state.jq), "FontSize", FS, "FontWeight","normal");
        grid(axM,"on");
        set(axM, "FontSize", FS, "LineWidth", 1.0, "TickDir","out", "Box","on");
    else
        % default behavior: plot vs mu and mark mu slider
        hold(axM,"on");
        okm = isfinite(mu);
        scatter(axM, mu(okm), chi(okm), 18, "filled", "MarkerFaceAlpha", 0.18, "MarkerEdgeAlpha", 0.18);
        [xm, ym_mean, ~, ~] = group_mean_std(mu(okm), chi(okm));
        if ~isempty(xm)
            plot(axM, xm, ym_mean, "-", "LineWidth", 1.8);
        end

        if isfinite(mu_pick) && ~isempty(xm)
            [~, jj] = min(abs(xm - mu_pick));
            plot(axM, xm(jj), ym_mean(jj), "o", "MarkerSize", 8, "LineWidth", 1.5);
        end
        hold(axM,"off");

        xlabel(axM, "\mu (eV)", "FontSize", FS);
        ylabel(axM, "Re \chi(iq,jq)", "FontSize", FS);
        title(axM, sprintf("Re \\chi vs \\mu @ T=%.6g, polar=%.6g (iq=%d,jq=%d)", ...
            Tnow, Pnow, state.iq, state.jq), "FontSize", FS, "FontWeight","normal");
        grid(axM,"on");
        set(axM, "FontSize", FS, "LineWidth", 1.0, "TickDir","out", "Box","on");
    end
end

end % ===== end main =====


% ============================================================
% helpers (ALL included)
% ============================================================

function pol = parse_polar_from_path(p)
    % Parse "polar" value from folder token like "polar_meV0.000"
    % Return NaN if not found.
    pol = NaN;
    s = char(p);
    tok = regexp(s, "polar[_\-]meV([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)", "tokens", "once");
    if ~isempty(tok)
        pol = str2double(tok{1});
    end
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

    % Case B: one extra run_tag layer (e.g. polar_meV0.000)
    sub = dir(D_dir);
    sub = sub([sub.isdir]);
    sub = sub(~ismember({sub.name},{'.','..'}));

    for i = 1:numel(sub)
        cand = fullfile(D_dir, sub(i).name);
        T2 = dir(fullfile(cand, "T*"));
        T2 = T2([T2.isdir]);
        T2 = T2(~ismember({T2.name},{'.','..'}));
        if ~isempty(T2)
            chi_run_dir = string(D_dir); % keep at D_dir, because we want multiple polars
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
    % Parse header lines starting with '#'
    % Accept:
    %   T =, T_K =
    %   doping =
    %   mu =, mu_eV =, EF =, EF_eV =
    %   mode = mu / doping / ...
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