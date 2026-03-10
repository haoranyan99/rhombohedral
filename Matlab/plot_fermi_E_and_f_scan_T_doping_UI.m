function out_png = plot_fermi_E_and_f_scan_T_doping_UI()
% plot_fermi_E_and_f_scan_T_doping_UI
% UI sliders for discrete (T, doping) points from debug fermiPatch outputs.
%
% Expected folder structure:
%   fermi_*/D*/T*/doping*/fermiPatch_*.txt
%
% Expected columns in each txt (debug mode):
%   idx iq jq kx ky occ_k_avg  E0..E{dim-1}  f0..f{dim-1}
%
% UI:
%   - Slider-1: T index
%   - Slider-2: doping index
%   - Edit box: band index n (0-based)
%   - 2-panel: En(k) heatmap & fn(k) heatmap
%   - Click on any pixel -> datatip shows idx, iq/jq, kx/ky, En, fn (for current band)
%
% Updates (vs your version):
%   - Use surf+view(2) instead of pcolor (better “cell center” alignment)
%   - Enable datacursor tooltip: idx + eigenvalue + occupancy at clicked point
%   - Top info text no longer shows file path

    out_png = "";

    % -----------------------------
    % choose root directory (e.g. fermi_2/D0.067)
    % -----------------------------
    start_dir = ".";
    root_dir = uigetdir(start_dir, "Select root directory like fermi_2/D0.067");
    if isequal(root_dir,0)
        error("No directory selected.");
    end
    root_dir = string(root_dir);

    % -----------------------------
    % scan files
    % -----------------------------
    DB = scan_fermi_debug_tree(root_dir);
    if isempty(DB) || isempty(DB.Tvals)
        error("No fermiPatch debug files found under: %s", root_dir);
    end

    % initial indices
    it = 1;
    id = 1;

    % -----------------------------
    % typography
    % -----------------------------
    FS = 14;
    set(groot, "defaultAxesTickLabelInterpreter", "latex");
    set(groot, "defaultLegendInterpreter", "latex");

    % -----------------------------
    % build figure
    % -----------------------------
    fig = figure("Color","w", "Units","pixels", "Position",[80 80 1360 720]);

    hInfo = uicontrol(fig, "Style","text", "Units","normalized", ...
        "Position",[0.05 0.965 0.93 0.030], ...
        "String","", "BackgroundColor","w", ...
        "FontSize", 11, "HorizontalAlignment","left");

    tl = tiledlayout(fig, 1, 2, "TileSpacing","compact", "Padding","compact");
    tl.Position = [0.06 0.18 0.92 0.75];
    ax1 = nexttile(tl, 1);
    ax2 = nexttile(tl, 2);

    % -----------------------------
    % sliders region
    % -----------------------------
    uicontrol(fig, "Style","text", "Units","normalized", ...
        "Position",[0.06 0.115 0.08 0.03], ...
        "String","T idx", "BackgroundColor","w", "FontSize", 11, ...
        "HorizontalAlignment","left");

    hT = uicontrol(fig, "Style","slider", "Units","normalized", ...
        "Position",[0.13 0.118 0.60 0.025], ...
        "Min",1, "Max",numel(DB.Tvals), "Value",it, ...
        "SliderStep", slider_steps(numel(DB.Tvals)), ...
        "Callback", @onTChanged);

    uicontrol(fig, "Style","text", "Units","normalized", ...
        "Position",[0.06 0.075 0.08 0.03], ...
        "String","dop idx", "BackgroundColor","w", "FontSize", 11, ...
        "HorizontalAlignment","left");

    hD = uicontrol(fig, "Style","slider", "Units","normalized", ...
        "Position",[0.13 0.078 0.60 0.025], ...
        "Min",1, "Max",numel(DB.byT{it}.xvals), "Value",id, ...
        "SliderStep", slider_steps(numel(DB.byT{it}.doping)), ...
        "Callback", @onDChanged);

    uicontrol(fig, "Style","text", "Units","normalized", ...
        "Position",[0.75 0.115 0.08 0.03], ...
        "String","band n", "BackgroundColor","w", "FontSize", 11, ...
        "HorizontalAlignment","left");

    hBand = uicontrol(fig, "Style","edit", "Units","normalized", ...
        "Position",[0.83 0.118 0.06 0.030], ...
        "String","0", "FontSize", 12, ...
        "Callback", @onBandChanged);

    uicontrol(fig, "Style","pushbutton", "Units","normalized", ...
        "Position",[0.91 0.114 0.07 0.038], ...
        "String","Apply", "FontSize", 12, ...
        "Callback", @onBandChanged);

    uicontrol(fig, "Style","pushbutton", "Units","normalized", ...
        "Position",[0.83 0.065 0.15 0.045], ...
        "String","Save PNG", "FontSize", 12, ...
        "Callback", @onSave);

    % -----------------------------
    % state
    % -----------------------------
    state = struct();
    state.root_dir = root_dir;
    state.DB = DB;
    state.it = it;
    state.id = id;
    state.n  = int32(0);
    state.in_path = "";
    state.meta = [];
    state.D = [];

    % datacursor: global (works for both axes)
    dcm = datacursormode(fig);
    set(dcm, "Enable","on", "DisplayStyle","window", "SnapToDataVertex","off");
    set(dcm, "UpdateFcn", @data_tip_update);

    % initial load & plot
    reload_current_file();
    update_plot();

    % =============================
    % callbacks
    % =============================
    function onTChanged(~,~)
        state.it = clamp_int(round(get(hT,"Value")), 1, numel(state.DB.Tvals));
        set(hT, "Value", state.it);

        ndop = numel(state.DB.byT{state.it}.doping);
        if ndop < 1, ndop = 1; end

        state.id = clamp_int(state.id, 1, ndop);
        set(hD, "Min",1, "Max",ndop, "Value",state.id, ...
            "SliderStep", slider_steps(ndop));

        reload_current_file();
        update_plot();
    end

    function onDChanged(~,~)
        ndop = numel(state.DB.byT{state.it}.doping);
        state.id = clamp_int(round(get(hD,"Value")), 1, max(1,ndop));
        set(hD, "Value", state.id);

        reload_current_file();
        update_plot();
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
        if strlength(state.in_path) == 0
            return;
        end
        [fdir, base, ~] = fileparts(state.in_path);
        out_png = fullfile(string(fdir), base + "_En_fn_Tdop_UI_n" + string(state.n) + ".png");
        exportgraphics(fig, out_png, "Resolution", 300);
        fprintf("Saved:\n  %s\n", out_png);
    end

    % =============================
    % core ops
    % =============================
    function reload_current_file()
        Tval = state.DB.Tvals(state.it);
        entryT = state.DB.byT{state.it};

        if isempty(entryT.files)
            state.in_path = "";
            state.meta = [];
            state.D = [];
            set(hInfo,"String", sprintf("No files under T=%gK", Tval));
            return;
        end

        ndop = numel(entryT.doping);
        if ndop < 1
            state.in_path = "";
            state.meta = [];
            state.D = [];
            set(hInfo,"String", sprintf("No doping folders under T=%gK", Tval));
            return;
        end
        state.id = clamp_int(state.id, 1, ndop);

        state.in_path = entryT.files(state.id);

        meta = parse_header_meta_fermi_debug_single(state.in_path);
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

        n = double(state.n) + 1; % MATLAB index
        EG = state.D.EG(:,:,n);
        fG = state.D.fG(:,:,n);

        maskE = isfinite(EG) & isfinite(state.D.kx0G) & isfinite(state.D.ky0G);
        maskF = isfinite(fG) & isfinite(state.D.kx0G) & isfinite(state.D.ky0G);

        dop0  = state.D.meta.doping;
        fill0 = state.D.meta.filling;
        EF0   = state.D.meta.EF_eV;
        T0    = state.D.meta.T_K;

        Emin = NaN; Emax = NaN;
        fmin = NaN; fmax = NaN;
        if any(maskE,'all')
            Emin = min(EG(maskE)); Emax = max(EG(maskE));
        end
        if any(maskF,'all')
            fmin = min(fG(maskF)); fmax = max(fG(maskF));
        end

        % NOTE: no file path here
        set(hInfo, "String", sprintf( ...
            "T=%.6g K | dop=%.6g (1e12 cm^-2) | fill=%.16g | EF=%.10g eV | band n=%d | En[min,max]=[%.6g,%.6g] | fn[min,max]=[%.6g,%.6g] | (click pixel to see idx/En/fn)", ...
            T0, dop0, fill0, EF0, state.n, Emin, Emax, fmin, fmax));

        % ---------- LEFT: En (use surf to align centers) ----------
        Z1 = EG; Z1(~maskE) = NaN;
        hS1 = surf(ax1, state.D.kx0G, state.D.ky0G, Z1, "EdgeColor","none");
        view(ax1, 2);
        axis(ax1,"equal"); axis(ax1,"tight");
        xlabel(ax1, '$(k_x-k_x^0)\;(1/\AA)$', "FontSize", FS);
        ylabel(ax1, '$(k_y-k_y^0)\;(1/\AA)$', "FontSize", FS);
        title(ax1, sprintf('$E_{%d}(\\mathbf{k})$ (eV)', state.n), ...
            "FontSize", FS, "FontWeight","normal");
        set(ax1, "FontSize", FS, "LineWidth", 1.0, "TickDir","out", "Box","on");
        colormap(ax1, jet(256));
        cb1 = colorbar(ax1, 'northoutside');
        cb1.Box = 'on';

        hold(ax1,"on");
        draw_BZ_rhombus_boundary(ax1, state.D);
        hold(ax1,"off");

        % store which panel this surface belongs to (for datatip)
        setappdata(hS1, "panel_tag", "E");

        % ---------- RIGHT: fn ----------
        Z2 = fG; Z2(~maskF) = NaN;
        hS2 = surf(ax2, state.D.kx0G, state.D.ky0G, Z2, "EdgeColor","none");
        view(ax2, 2);
        axis(ax2,"equal"); axis(ax2,"tight");
        xlabel(ax2, '$(k_x-k_x^0)\;(1/\AA)$', "FontSize", FS);
        ylabel(ax2, '$(k_y-k_y^0)\;(1/\AA)$', "FontSize", FS);
        title(ax2, sprintf('$f_{%d}(\\mathbf{k})$', state.n), ...
            "FontSize", FS, "FontWeight","normal");
        set(ax2, "FontSize", FS, "LineWidth", 1.0, "TickDir","out", "Box","on");
        colormap(ax2, flipud(hot(256)));   % big value -> dark
        caxis(ax2, [0, 1]);
        cb2 = colorbar(ax2, 'northoutside');
        cb2.Box = 'on';

        hold(ax2,"on");
        draw_BZ_rhombus_boundary(ax2, state.D);
        hold(ax2,"off");

        setappdata(hS2, "panel_tag", "f");
    end

    % =============================
    % datatip update function
    % =============================
    function txt = data_tip_update(~, event)
        % event.Position is [x y z] in plotted coords (kx0, ky0, value)
        if isempty(state.D) || ~state.D.has_debug_bands
            txt = {"(no data)"};
            return;
        end

        pos = event.Position;
        xq = pos(1);
        yq = pos(2);

        % find nearest grid point (robust even if there are NaNs)
        X = state.D.kx0G;
        Y = state.D.ky0G;

        good = isfinite(X) & isfinite(Y);
        if ~any(good,'all')
            txt = {"(no valid grid)"};
            return;
        end

        dx = X - xq;
        dy = Y - yq;
        dist2 = dx.^2 + dy.^2;
        dist2(~good) = +inf;
        [~, lin] = min(dist2(:));
        [rr, cc] = ind2sub(size(dist2), lin);

        idx0 = NaN;
        if isfield(state.D, "idxG") && ~isempty(state.D.idxG)
            idx0 = state.D.idxG(rr,cc);
        end

        iq0 = NaN; jq0 = NaN;
        if isfield(state.D,"iqG"), iq0 = state.D.iqG(rr,cc); end
        if isfield(state.D,"jqG"), jq0 = state.D.jqG(rr,cc); end

        kx0 = X(rr,cc);
        ky0 = Y(rr,cc);

        n = double(state.n) + 1;
        En = state.D.EG(rr,cc,n);
        fn = state.D.fG(rr,cc,n);

        % panel tag (E or f), optional
        tag = "";
        try
            h = event.Target;
            tag = string(getappdata(h, "panel_tag"));
        catch
        end

        % datatip strings (NO path)
        txt = {
            sprintf("idx = %g", idx0)
            sprintf("(iq,jq) = (%g, %g)", iq0, jq0)
            sprintf("(kx,ky) = (%.8g, %.8g)", kx0, ky0)
            sprintf("band n = %d", state.n)
            sprintf("E_n = %.10g eV", En)
            sprintf("f_n = %.10g", fn)
        };

        if strlength(tag) > 0
            txt{end+1} = sprintf("panel = %s", tag);
        end
    end
end

% =========================================================
% helpers: scanning
% =========================================================
function DB = scan_fermi_debug_tree(root_dir)
    root_dir = char(root_dir);
    if ~exist(root_dir,"dir")
        error("scan: root_dir not found: %s", root_dir);
    end

    files = dir(fullfile(root_dir, "**", "fermiPatch*.txt"));
    if isempty(files)
        DB = struct("Tvals",[],"byT",{{}});
        return;
    end

    recs = struct("T",{}, "x",{}, "path",{}, "mode",{});
    Tlist = [];

    for i = 1:numel(files)
        fpath = fullfile(files(i).folder, files(i).name);

        % ---- read header (cheap) ----
        try
            meta = parse_header_meta_fermi_debug_single(fpath);
        catch
            continue;
        end

        % ---- T: prefer header, fallback to path /Txx/ ----
        Tval = meta.T_K;
        if ~isfinite(Tval)
            fp = strrep(fpath, "\", "/");
            tokT = regexp(fp, "/T([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)/", "tokens", "once");
            if ~isempty(tokT), Tval = str2double(tokT{1}); end
        end
        if ~isfinite(Tval), continue; end

        mode = string(meta.mode);
        % ---- x: slider-2 key (mu-mode -> EF, otherwise -> doping; fallback EF; fallback index) ----
        xval = NaN;
        if strcmpi(strtrim(mode), "mu")
            xval = meta.EF_eV;
        else
            xval = meta.doping;
            if ~isfinite(xval), xval = meta.EF_eV; end
        end
        if ~isfinite(xval)
            xval = i; % last resort: stable order
        end

        recs(end+1).T    = Tval; %#ok<AGROW>
        recs(end).x      = xval;
        recs(end).path   = string(fpath);
        recs(end).mode   = mode;

        Tlist(end+1) = Tval; %#ok<AGROW>
    end

    if isempty(recs)
        DB = struct("Tvals",[],"byT",{{}});
        return;
    end

    Tvals = unique(Tlist);
    Tvals = sort(Tvals);

    byT = cell(numel(Tvals),1);
    for it = 1:numel(Tvals)
        Tv = Tvals(it);
        mask = arrayfun(@(r) r.T==Tv, recs);
        rs = recs(mask);

        xs = [rs.x];
        [xs_sorted, ord] = sort(xs);
        rs = rs(ord);

        % store
        byT{it} = struct();
        byT{it}.T     = Tv;
        byT{it}.xvals = xs_sorted(:);               % <<<< IMPORTANT: not "doping" anymore
        byT{it}.files = string({rs.path}).';
        % if mixed modes exist, keep a tag for label (optional)
        byT{it}.mode  = (numel(unique(string({rs.mode})))==1) * string(rs(1).mode) + ...
                        (numel(unique(string({rs.mode})))>1)  * "mixed";
    end

    DB = struct();
    DB.Tvals = Tvals(:);
    DB.byT   = byT;
end


function st = slider_steps(N)
    if N <= 1
        st = [1 1];
    else
        st = [1/(N-1), min(10/(N-1),1)];
    end
end

function x = clamp_int(x, lo, hi)
    x = max(lo, min(hi, x));
end

% =========================================================
% helpers: file parsing & gridding
% =========================================================
function meta = parse_header_meta_fermi_debug_single(in_path)
    meta = struct();
    meta.Nk      = NaN;
    meta.dk_frac = NaN;
    meta.filling = NaN;
    meta.doping  = NaN;
    meta.T_K     = NaN;
    meta.Dfield  = NaN;
    meta.EF_eV   = NaN;

    meta.mode    = "";   % <<< NEW
    meta.model   = "";   % optional
    meta.b1 = [NaN; NaN];
    meta.b2 = [NaN; NaN];

    fid = fopen(in_path,'r');
    if fid < 0, error("Cannot open file: %s", in_path); end

    num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end
        s0 = strtrim(tline);
        if ~startsWith(s0, "#"), break; end

        % ---- b1/b2 vector lines ----
        tb1 = regexp(s0, "#\s*b1\s*=\s*\(\s*"+num+"\s*,\s*"+num+"\s*\)", "tokens", "once");
        if ~isempty(tb1)
            meta.b1 = [str2double(tb1{1}); str2double(tb1{2})];
        end
        tb2 = regexp(s0, "#\s*b2\s*=\s*\(\s*"+num+"\s*,\s*"+num+"\s*\)", "tokens", "once");
        if ~isempty(tb2)
            meta.b2 = [str2double(tb2{1}); str2double(tb2{2})];
        end

        % ---- generic "key = val" ----
        s = strtrim(erase(s0, "#"));
        tok = regexp(s, "^([A-Za-z0-9_\-\.\[\]\^]+)\s*=\s*(.*)$", "tokens", "once");
        if ~isempty(tok)
            key = strtrim(tok{1});
            val = strtrim(tok{2});

            % strip possible quotes
            if startsWith(val,'"') && endsWith(val,'"')
                val = extractBetween(string(val), 2, strlength(string(val))-1);
                val = char(val);
            end

            % numeric?
            tv = regexp(val, "^"+num+"$", "tokens", "once");
            if ~isempty(tv)
                vnum = str2double(tv{1});
                switch key
                    case {"filling","filling_target","fill"}
                        meta.filling = vnum;
                    case {"doping_1e12cm^-2","doping_1e12cm-2","doping"}
                        meta.doping = vnum;
                    case {"T_K","T"}
                        meta.T_K = vnum;
                    case {"Dfield_eV","Dfield"}
                        meta.Dfield = vnum;
                    case {"EF_eV","EF","mu","mu_eV"}
                        meta.EF_eV = vnum;
                end
            else
                % string keys
                switch key
                    case {"mode"}
                        meta.mode = string(val);
                    case {"model"}
                        meta.model = string(val);
                end
            end
        end

        % ---- Nk / dk_frac in kmesh line ----
        if contains(s, "Nk=") && ~isfinite(meta.Nk)
            tk = regexp(s, "Nk\s*=\s*(\d+)", "tokens", "once");
            if ~isempty(tk), meta.Nk = str2double(tk{1}); end
        end
        if contains(s, "dk_frac=") && ~isfinite(meta.dk_frac)
            td = regexp(s, "dk_frac\s*=\s*"+num, "tokens", "once");
            if ~isempty(td), meta.dk_frac = str2double(td{1}); end
        end
    end
    fclose(fid);
end


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

    % origin shift (kx0,ky0)
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
