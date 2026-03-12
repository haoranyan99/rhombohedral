function plot_chi_vs_mu()
% plot_chi_vs_mu_UI
% UI sliders to browse chiEf outputs and plot Re(chi) vs mu at fixed (T, iq, jq).
%
% Expected folder structure (recommended):
%   data/chiEf_{suffix}/D{Dfield}/T{T_K}/mu{mu}/chiEf_*.txt
%
% The script is tolerant to an extra run_tag layer:
%   D_dir/T*
%   or D_dir/*/T*
%
% It reads ONE chi file to infer available iq/jq lists, then plots
% Re(chi) at selected (iq,jq) vs mu across mu folders.

    % -----------------------------
    % choose D directory
    % -----------------------------
    base_dir = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\chi_sk_mu_500\D0.067";
    if ~exist(base_dir,"dir"), base_dir = "."; end

    D_dir = uigetdir(base_dir, "Select ONE D folder (e.g. data/chiEf_6L_3/D0.000)");
    if D_dir == 0
        fprintf("Canceled.\n");
        return;
    end
    D_dir = string(D_dir);
    fprintf("Selected D folder:\n  %s\n", D_dir);

    % -----------------------------
    % resolve where T folders live
    %   prefer: D_dir/T*
    %   else:   D_dir/*/T*
    % -----------------------------
    [chi_run_dir, Tlist] = resolve_T_root(D_dir);
    fprintf("Resolved T root:\n  %s\n", chi_run_dir);

    % parse T values
    T_vals = nan(numel(Tlist),1);
    for i = 1:numel(Tlist)
        T_vals(i) = sscanf(Tlist(i).name, "T%f");
    end
    [T_vals, tidx] = sort(T_vals);
    Tlist = Tlist(tidx);

    % -----------------------------
    % infer mu list from FIRST T folder
    % -----------------------------
    mu_dirs = dir(fullfile(chi_run_dir, Tlist(1).name, "mu*"));
    mu_dirs = mu_dirs([mu_dirs.isdir]);
    if isempty(mu_dirs)
        error("No mu folders under: %s", fullfile(chi_run_dir, Tlist(1).name));
    end

    mu_vals = nan(numel(mu_dirs),1);
    for i = 1:numel(mu_dirs)
        mu_vals(i) = sscanf(mu_dirs(i).name, "mu%f");
    end
    [mu_vals, midx] = sort(mu_vals);
    mu_dirs = mu_dirs(midx);

    % -----------------------------
    % read one chi file to infer iq/jq lists
    % -----------------------------
    chi0 = find_latest_chi_in_folder(fullfile(chi_run_dir, Tlist(1).name, mu_dirs(1).name));
    if strlength(chi0)==0
        error("No chi_*.txt (or chi*.txt) found in: %s", ...
              fullfile(chi_run_dir, Tlist(1).name, mu_dirs(1).name));
    end
    D0 = read_chi_file_auto(chi0);

    iq_list = sort(unique(D0.iq));
    jq_list = sort(unique(D0.jq));

    % defaults
    T_idx  = 1;
    iq_idx = pick_default_index(iq_list, 0);
    jq_idx = pick_default_index(jq_list, 0);

    % -----------------------------
    % UI
    % -----------------------------
    FS = 16;
    fig = figure("Color","w","Units","pixels","Position",[80 80 1200 820]);

    ax = axes(fig,"Position",[0.10 0.32 0.85 0.62]);
    box(ax,"on"); grid(ax,"on"); hold(ax,"on");
    set(ax,"FontSize",FS,"TickDir","out","LineWidth",1.0);

    hLine = plot(ax, nan, nan, "o-", "LineWidth",1.8,"MarkerSize",6);

    xlabel(ax,"\mu (eV)","FontSize",FS);
    ylabel(ax,"Re(\chi) (arb.)","FontSize",FS);

    % ---- sliders ----
    sT = uicontrol(fig,"Style","slider","Units","normalized", ...
        "Position",[0.10 0.22 0.85 0.04], ...
        "Min",1,"Max",numel(T_vals),"Value",T_idx, ...
        "SliderStep",slider_step(numel(T_vals)));

    sIQ = uicontrol(fig,"Style","slider","Units","normalized", ...
        "Position",[0.10 0.14 0.85 0.04], ...
        "Min",1,"Max",numel(iq_list),"Value",iq_idx, ...
        "SliderStep",slider_step(numel(iq_list)));

    sJQ = uicontrol(fig,"Style","slider","Units","normalized", ...
        "Position",[0.10 0.06 0.85 0.04], ...
        "Min",1,"Max",numel(jq_list),"Value",jq_idx, ...
        "SliderStep",slider_step(numel(jq_list)));

    mkLabel(fig, FS, "T",  0.22);
    mkLabel(fig, FS, "iq", 0.14);
    mkLabel(fig, FS, "jq", 0.06);

    txtT  = mkValue(fig, FS, 0.26);
    txtIQ = mkValue(fig, FS, 0.18);
    txtJQ = mkValue(fig, FS, 0.10);

    set([sT,sIQ,sJQ],"Callback",@onSlide);

    update_plot();

    % =============================
    function onSlide(~,~)
        sT.Value  = round(sT.Value);
        sIQ.Value = round(sIQ.Value);
        sJQ.Value = round(sJQ.Value);
        update_plot();
    end

    function update_plot()
        it  = clampi(round(sT.Value),  1, numel(T_vals));
        ii  = clampi(round(sIQ.Value), 1, numel(iq_list));
        jj  = clampi(round(sJQ.Value), 1, numel(jq_list));

        T_now  = T_vals(it);
        iq_now = iq_list(ii);
        jq_now = jq_list(jj);

        txtT.String  = sprintf("T = %.4f K", T_now);
        txtIQ.String = sprintf("iq = %d", iq_now);
        txtJQ.String = sprintf("jq = %d", jq_now);

        re_vec = nan(numel(mu_vals),1);
        qx_now = NaN; qy_now = NaN;

        T_folder = fullfile(chi_run_dir, Tlist(it).name);

        for k = 1:numel(mu_vals)
            mu_folder = fullfile(T_folder, sprintf("mu%.6f", mu_vals(k)));

            % tolerate tiny format mismatch in folder name
            if ~exist(mu_folder,"dir")
                mu_folder = find_matching_mu_folder(T_folder, mu_vals(k));
                if strlength(mu_folder)==0, continue; end
            end

            chi_file = find_latest_chi_in_folder(mu_folder);
            if strlength(chi_file)==0, continue; end

            D = read_chi_file_auto(chi_file);
            m = (D.iq==iq_now) & (D.jq==jq_now);
            if any(m)
                id = find(m,1,"first");
                re_vec(k) = real(D.chi(id));
                if ~isfinite(qx_now)
                    qx_now = D.qx(id);
                    qy_now = D.qy(id);
                end
            end
        end

        ok = isfinite(re_vec);
        set(hLine,"XData",mu_vals(ok),"YData",re_vec(ok));

        if isfinite(qx_now)
            title(ax, sprintf( ...
                "Re\\chi vs \\mu | T=%.4f K | (iq,jq)=(%d,%d) | (qx,qy)=(%.6g, %.6g)", ...
                T_now, iq_now, jq_now, qx_now, qy_now), ...
                "Interpreter","tex","FontWeight","normal");
        else
            title(ax, sprintf( ...
                "Re\\chi vs \\mu | T=%.4f K | (iq,jq)=(%d,%d)", ...
                T_now, iq_now, jq_now), ...
                "Interpreter","tex","FontWeight","normal");
        end
    end
end

% ============================================================
% Resolve where T folders live:
%   prefer: D_dir/T*
%   else:   D_dir/*/T*
% ============================================================
function [chi_run_dir, Tlist] = resolve_T_root(D_dir)
    % case A
    Tlist = dir(fullfile(D_dir, "T*"));
    Tlist = Tlist([Tlist.isdir]);
    if ~isempty(Tlist)
        chi_run_dir = D_dir;
        return;
    end

    % case B: search one level down
    sub = dir(D_dir);
    sub = sub([sub.isdir]);
    sub = sub(~ismember({sub.name},{'.','..'}));

    for i = 1:numel(sub)
        cand = fullfile(D_dir, sub(i).name);
        T2 = dir(fullfile(cand, "T*"));
        T2 = T2([T2.isdir]);
        if ~isempty(T2)
            chi_run_dir = string(cand);
            Tlist = T2;
            return;
        end
    end

    error("Cannot find any T* folders under:\n  %s\nor under one subfolder level.", D_dir);
end

function mu_folder = find_matching_mu_folder(T_folder, mu_val)
    mu_folder = "";
    dd = dir(fullfile(T_folder, "mu*"));
    dd = dd([dd.isdir]);
    for k = 1:numel(dd)
        v = sscanf(dd(k).name, "mu%f");
        if ~isempty(v) && abs(v - mu_val) < 5e-7  % mu uses 6 decimals typically
            mu_folder = string(fullfile(T_folder, dd(k).name));
            return;
        end
    end
end

function p = find_latest_chi_in_folder(folder)
% Prefer chiEf_*.txt; fallback chi_*.txt
    p = "";
    if ~exist(folder,"dir"), return; end

    files = dir(fullfile(folder,"chi_*.txt"));
    if isempty(files)
        files = dir(fullfile(folder,"chi*.txt"));
    end
    if isempty(files), return; end

    [~,ord] = sort({files.name});
    p = string(fullfile(folder,files(ord(end)).name));
end

function D = read_chi_file_auto(fname)
% Supports numeric formats:
%   (A) 7 cols:  iq jq qx qy Re Im nKpair
%   (B) 8 cols:  idx iq jq qx qy Re Im nKpair
%   (C) 9 cols:  idx iq jq qx qy Re Im nKpair nNesting
%
% If nNesting exists, D.nNest is filled, otherwise D.nNest = NaN.

    fid = fopen(fname, 'r');
    if fid < 0, error("Cannot open %s", fname); end

    % skip header
    while true
        pos = ftell(fid);
        line = fgetl(fid);
        if ~ischar(line), break; end
        t = strtrim(line);
        if isempty(t), continue; end
        if startsWith(t,"#")
            continue;
        else
            fseek(fid, pos, 'bof');
            break;
        end
    end

    pos0 = ftell(fid);
    first = fgetl(fid);
    if ~ischar(first)
        fclose(fid);
        error("No numeric data in %s", fname);
    end
    nums = sscanf(first, '%f');
    ncol = numel(nums);
    fseek(fid, pos0, 'bof');

    D.nNest = [];

    if ncol == 9
        C = textscan(fid, '%f %f %f %f %f %f %f %f %f', 'CollectOutput', true);
        M = C{1};
        D.iq    = M(:,2);
        D.jq    = M(:,3);
        D.qx    = M(:,4);
        D.qy    = M(:,5);
        D.chi   = M(:,6) + 1i*M(:,7);
        D.nK    = M(:,8);
        D.nNest = M(:,9);
    elseif ncol == 8
        C = textscan(fid, '%f %f %f %f %f %f %f %f', 'CollectOutput', true);
        M = C{1};
        D.iq  = M(:,2);
        D.jq  = M(:,3);
        D.qx  = M(:,4);
        D.qy  = M(:,5);
        D.chi = M(:,6) + 1i*M(:,7);
        D.nK  = M(:,8);
        D.nNest = nan(size(D.nK));
    elseif ncol == 7
        C = textscan(fid, '%f %f %f %f %f %f %f', 'CollectOutput', true);
        M = C{1};
        D.iq  = M(:,1);
        D.jq  = M(:,2);
        D.qx  = M(:,3);
        D.qy  = M(:,4);
        D.chi = M(:,5) + 1i*M(:,6);
        D.nK  = M(:,7);
        D.nNest = nan(size(D.nK));
    else
        fclose(fid);
        error("Unexpected numeric column count (%d) in %s", ncol, fname);
    end

    fclose(fid);
end

function idx = pick_default_index(list,target)
    k = find(list==target,1);
    if isempty(k), idx = round((numel(list)+1)/2);
    else, idx = k;
    end
end

function step = slider_step(n)
    if n<=1, step=[1 1];
    else, step=[1/(n-1) 5/(n-1)]; step(step>1)=1;
    end
end

function x = clampi(x,a,b)
    x = min(max(x,a),b);
end

function mkLabel(fig, FS, str, y)
    uicontrol(fig,"Style","text","Units","normalized", ...
        "Position",[0.02 y 0.07 0.04], ...
        "String",str,"FontSize",FS,"BackgroundColor","w");
end

function h = mkValue(fig, FS, y)
    h = uicontrol(fig,"Style","text","Units","normalized", ...
        "Position",[0.10 y 0.50 0.03], ...
        "String","","FontSize",FS, ...
        "BackgroundColor","w","HorizontalAlignment","left");
end
