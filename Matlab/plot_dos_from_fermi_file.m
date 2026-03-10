function out_png_list = plot_dos_from_fermi_file()
% plot_dos_from_fermi_file
% Select ONE debug fermiPatch_*.txt and compute DOS from all eigenvalues E0..E11.
% DOS is computed ONLY within a user-defined energy window [Emin_show, Emax_show],
% using bin width dE = (Emax_show - Emin_show) / Nbin.

    out_png_list = strings(0);

    % ============================================================
    % USER: DISPLAY RANGE + GRID
    % ============================================================
    range = struct();

    % choose ONE mode:
    range.mode = "absolute";      % "EF_window" or "absolute"

    % (A) EF_window: [EF - dEwin, EF + dEwin]
    range.dEwin = 0.05;            % eV

    % (B) absolute: [Emin_show, Emax_show]
    range.Emin_show = 0.90;        % eV
    range.Emax_show = 1.1;        % eV

    range.clip_to_data = true;     % clip window into [min(Eall), max(Eall)]
    range.Nbin = 1200;             % <<< N in your "dE = range/N"

    % ============================================================
    % USER: BROADENING
    % ============================================================
    par = struct();
    par.eta_eV = 0.002;            % Gaussian sigma (eV)

    % -----------------------------
    % choose file
    % -----------------------------
    default_path = fullfile("data", "*.txt");
    [fname, fpath] = uigetfile({'*.txt','fermiPatch (*.txt)'}, ...
                               'Select fermiPatch_*.txt file', default_path);
    if isequal(fname,0), error("No file selected."); end
    in_path = string(fullfile(fpath, fname));
    if exist(char(in_path),'file') ~= 2
        error("Input file not found: %s", in_path);
    end
    [~, base, ~] = fileparts(in_path);

    % -----------------------------
    % parse header meta
    % -----------------------------
    meta = parse_header_meta_fermiPatch_safe(in_path);

    % -----------------------------
    % read numeric block robustly
    % -----------------------------
    [M, ncol] = read_numeric_block_textscan_safe(in_path, meta.columns_line);

    % -----------------------------
    % extract eigenvalues E0..E{Nb-1}
    % -----------------------------
    nBase = 6; % idx iq jq kx ky occ_k_avg
    if ncol <= nBase
        error("No eigenvalue columns found. ncol=%d in %s", ncol, in_path);
    end

    rest = ncol - nBase;
    if mod(rest,2)==0
        Nb = rest/2;     % debug format: E(Nb)+f(Nb)
    else
        Nb = rest;       % non-debug: only E
    end

    Ecols = (nBase+1):(nBase+Nb);
    Eall = M(:, Ecols);
    Eall = Eall(:);
    Eall = Eall(isfinite(Eall));
    if isempty(Eall)
        error("No valid eigenvalues extracted from: %s", in_path);
    end

    Emin_data = min(Eall);
    Emax_data = max(Eall);

    % -----------------------------
    % choose [Emin_show, Emax_show]
    % -----------------------------
    if range.mode == "EF_window"
        if ~isfinite(meta.EF)
            error("range.mode='EF_window' but EF not found in header: %s", in_path);
        end
        Emin_show = meta.EF - range.dEwin;
        Emax_show = meta.EF + range.dEwin;
    elseif range.mode == "absolute"
        Emin_show = range.Emin_show;
        Emax_show = range.Emax_show;
    else
        error("Unknown range.mode: %s", range.mode);
    end

    if range.clip_to_data
        Emin_show = max(Emin_show, Emin_data);
        Emax_show = min(Emax_show, Emax_data);
    end
    if ~(Emax_show > Emin_show)
        error("Invalid display range: [%.6g, %.6g]", Emin_show, Emax_show);
    end

    % ============================================================
    % ONLY compute DOS within [Emin_show, Emax_show]
    % dE = (range)/Nbin  (your requirement)
    % ============================================================
    Nbin = range.Nbin;
    if ~(isscalar(Nbin) && Nbin>=10 && isfinite(Nbin))
        error("range.Nbin must be a finite scalar >= 10");
    end

    width = (Emax_show - Emin_show);
    dE = width / Nbin;

    edges = linspace(Emin_show, Emax_show, Nbin+1);
    Ecent = Emin_show + ( (0.5:1:(Nbin-0.5))' ) * dE;   % bin centers

    % keep only eigenvalues inside window
    inwin = (Eall >= Emin_show) & (Eall <= Emax_show);
    Ewin = Eall(inwin);

    if isempty(Ewin)
        error("No eigenvalues lie within the chosen window [%.6g, %.6g].", Emin_show, Emax_show);
    end

    counts = histcounts(Ewin, edges);
    dos_raw = (counts(:) / numel(Eall)) / dE;
    % NOTE: normalize by total states numel(Eall) (so DOS is comparable across windows)
    % If you prefer DOS normalized within window only, replace numel(Eall) -> numel(Ewin).

    % Gaussian smoothing on this window grid (edge effects possible)
    sigma = par.eta_eV;
    L = max(3, ceil(4*sigma/dE));
    x = (-L:L)' * dE;
    ker = exp(-(x.^2)/(2*sigma^2));
    ker = ker / sum(ker);

    DOS = conv(dos_raw, ker, 'same');

    % -----------------------------
    % plot
    % -----------------------------
    out_dir = "plot";
    if ~exist(out_dir,"dir"), mkdir(out_dir); end

    FS = 14;
    fig = figure('Color','w', 'Units','pixels', 'Position',[120 120 640 440]);
    ax = axes(fig); hold(ax,'on');

    plot(ax, Ecent, DOS, 'LineWidth', 2.0);
    xlabel(ax, 'E (eV)', 'Interpreter','none', 'FontSize',FS);
    ylabel(ax, 'DOS (arb., normalized)', 'Interpreter','none', 'FontSize',FS);

    ttl = make_title_safe(meta, par, sprintf('DOS in [%.4f, %.4f] eV, dE=%.3g, Nb=%d', ...
        Emin_show, Emax_show, dE, Nb));
    title(ax, ttl, 'Interpreter','none', 'FontSize',FS, 'FontWeight','normal');

    set(ax, 'FontSize',FS, 'LineWidth',1.0, 'TickDir','out', 'Box','on');
    grid(ax,'on');
    xlim(ax, [Emin_show Emax_show]);

    if isfinite(meta.EF) && meta.EF >= Emin_show && meta.EF <= Emax_show
        xline(ax, meta.EF, '--', 'LineWidth', 1.0);
    end

    out_png = fullfile(out_dir, base + "_DOS_vsE.png");
    exportgraphics(fig, out_png, 'Resolution',300);
    fprintf("Saved figure:\n  %s\n", out_png);
    out_png_list(end+1) = string(out_png);
end

% ============================================================
% Safe header parser (no regexp directly used as logical)
% ============================================================
function meta = parse_header_meta_fermiPatch_safe(in_path)
    meta = struct();
    meta.model      = "";
    meta.nLayer     = NaN;
    meta.Dfield_eV  = NaN;
    meta.T_K        = NaN;
    meta.doping     = NaN;
    meta.EF         = NaN;
    meta.kmesh_type = "";
    meta.Nk         = NaN;
    meta.dk_frac    = NaN;
    meta.columns_line = "";

    fid = fopen(char(in_path),'r');
    if fid < 0, error("Cannot open file: %s", in_path); end
    c = onCleanup(@() fclose(fid));

    num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end
        s0 = strtrim(tline);
        if ~startsWith(s0, "#"), break; end

        s = strtrim(erase(s0, "#"));

        if strlength(meta.columns_line)==0 && startsWith(lower(strtrim(s)), "columns:")
            meta.columns_line = string(s);
        end

        if strlength(meta.model)==0 && contains(s, "model")
            tok = regexp(s, "(^|\s)model\s*=\s*([A-Za-z0-9_]+)", "tokens", "once");
            if ~isempty(tok), meta.model = string(tok{2}); end
        end

        if ~isfinite(meta.nLayer) && contains(s, "nLayer")
            tok = regexp(s, "(^|\s)nLayer\s*=\s*(\d+)", "tokens", "once");
            if ~isempty(tok), meta.nLayer = str2double(tok{2}); end
        end

        if ~isfinite(meta.Dfield_eV) && contains(s, "Dfield_eV")
            tok = regexp(s, "(^|\s)Dfield_eV\s*=\s*" + num, "tokens", "once");
            if ~isempty(tok), meta.Dfield_eV = str2double(tok{2}); end
        end

        if ~isfinite(meta.Dfield_eV) && contains(s, "Dfield")
            tok = regexp(s, "(^|\s)Dfield\s*=\s*" + num, "tokens", "once");
            if ~isempty(tok), meta.Dfield_eV = str2double(tok{2}); end
        end

        if ~isfinite(meta.T_K) && ~isempty(regexp(s, "^\s*T\s*=", "once"))
            tok = regexp(s, "T\s*=\s*" + num, "tokens", "once");
            if ~isempty(tok), meta.T_K = str2double(tok{1}); end
        end

        if ~isfinite(meta.doping) && ~isempty(regexp(s, "(^|\s)doping\s*=", "once"))
            tok = regexp(s, "(^|\s)doping\s*=\s*" + num, "tokens", "once");
            if ~isempty(tok), meta.doping = str2double(tok{2}); end
        end

        if ~isfinite(meta.EF) && ~isempty(regexp(s, "(^|\s)EF\s*=", "once"))
            tok = regexp(s, "(^|\s)EF\s*=\s*" + num, "tokens", "once");
            if ~isempty(tok), meta.EF = str2double(tok{2}); end
        end

        if (strlength(meta.kmesh_type)==0 || ~isfinite(meta.Nk) || ~isfinite(meta.dk_frac)) && contains(s, "kmesh")
            tok = regexp(s, "kmesh\s*=\s*([A-Za-z0-9_]+)", "tokens", "once");
            if ~isempty(tok), meta.kmesh_type = string(tok{1}); end

            tk = regexp(s, "Nk\s*=\s*(\d+)", "tokens", "once");
            if ~isempty(tk), meta.Nk = str2double(tk{1}); end

            td = regexp(s, "dk_frac\s*=\s*" + num, "tokens", "once");
            if ~isempty(td), meta.dk_frac = str2double(td{1}); end
        end
    end
end

% ============================================================
% Robust numeric reader via textscan('%f'), ignoring '#'
% ============================================================
function [M, ncol] = read_numeric_block_textscan_safe(in_path, columns_line)
    fid = fopen(char(in_path), 'r');
    if fid < 0, error("Cannot open: %s", in_path); end
    c = onCleanup(@() fclose(fid));

    % skip header comments
    pos = ftell(fid);
    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end
        s = strtrim(tline);
        if startsWith(s, "#") || strlength(s)==0
            pos = ftell(fid);
            continue;
        end
        fseek(fid, pos, 'bof');
        break;
    end

    data = textscan(fid, '%f', 'Delimiter', {' ','\t'}, 'MultipleDelimsAsOne', true);
    v = data{1};
    if isempty(v), error("No numeric data parsed from: %s", in_path); end

    % infer ncol from columns_line
    ncol = NaN;
    if strlength(columns_line) > 0
        s = char(columns_line);
        s = regexprep(s, "^[Cc]olumns:\s*", "");
        s = regexprep(s, "\s+", " ");
        toks = strsplit(strtrim(s), " ");
        if numel(toks) >= 7
            ncol = numel(toks);
        end
    end

    if isnan(ncol)
        cand = [30 28 26 24 22 20 18 16 14 12 10 8 7];
        for cc = cand
            if mod(numel(v), cc)==0
                ncol = cc; break;
            end
        end
    end
    if isnan(ncol)
        error("Cannot infer column count. Numeric token count=%d", numel(v));
    end

    nrow = floor(numel(v)/ncol);
    v = v(1:nrow*ncol);
    M = reshape(v, [ncol, nrow]).';
end

function ttl = make_title_safe(meta, par, line1)
    lines = {};
    if strlength(meta.model) > 0
        lines{end+1} = sprintf('%s  |  %s', upper(char(meta.model)), line1);
    else
        lines{end+1} = line1;
    end

    meta2 = {};
    if isfinite(meta.nLayer),    meta2{end+1} = sprintf('Nl=%d', meta.nLayer); end %#ok<AGROW>
    if isfinite(meta.Dfield_eV), meta2{end+1} = sprintf('D=%.6g eV', meta.Dfield_eV); end %#ok<AGROW>
    if isfinite(meta.T_K),       meta2{end+1} = sprintf('T=%.6g K', meta.T_K); end %#ok<AGROW>
    if isfinite(meta.doping),    meta2{end+1} = sprintf('doping=%.6g', meta.doping); end %#ok<AGROW>
    if isfinite(meta.EF),        meta2{end+1} = sprintf('EF=%.6g eV', meta.EF); end %#ok<AGROW>
    meta2{end+1} = sprintf('sigma=%.4g eV', par.eta_eV); %#ok<AGROW>

    lines{end+1} = strjoin(meta2, ', ');
    ttl = strjoin(lines, newline);
end
