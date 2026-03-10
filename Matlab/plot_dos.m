function out_png_list = plot_dos()
% plot_test_RG_dos
% Plot RG DOS from text file (NEW format: header carries metadata).
%
% Data block examples (comment lines start with '#'):
%   (A) i  E(eV)  DOS
%   (B) i  E(eV)  filling  DOS
%   (C) i  E(eV)  doping_1e12cm^-2  DOS
%   (D) i  E(eV)  filling  doping_1e12cm^-2  DOS
%
% This script makes up to 3 figures:
%   1) DOS vs E
%   2) DOS vs filling   (if exists)
%   3) DOS vs doping    (if exists)
%
% Output PNGs saved into ./plot with suffixes:
%   _vsE, _vsFilling, _vsDoping
%
% NOTE:
% - We DO NOT parse metadata from filename anymore.
% - We parse metadata (model, Nl, Dfield_eV, T_K, kmesh...) from header.

    out_png_list = strings(0);

    % -----------------------------
    % choose file
    % -----------------------------
    default_path = fullfile("data", "*.txt");
    [fname, fpath] = uigetfile({'*.txt','DOS data (*.txt)'}, ...
                               'Select RG DOS file', default_path);
    if isequal(fname,0)
        error("No file selected.");
    end
    in_path = fullfile(fpath, fname);
    if ~isfile(in_path)
        error("Input file not found: %s", in_path);
    end

    [~, base, ~] = fileparts(in_path);

    % -----------------------------
    % parse header meta + columns
    % -----------------------------
    meta = parse_header_meta_dos(in_path);

    % -----------------------------
    % read numeric block (skip comments)
    % -----------------------------
    raw = readmatrix(in_path, "FileType","text", "CommentStyle","#");
    if isempty(raw) || size(raw,2) < 3
        error("Bad DOS format: need at least 3 numeric columns.\nFile: %s", in_path);
    end

    % -----------------------------
    % decide columns by header first (most robust)
    % -----------------------------
    [E, filling, doping, DOS] = extract_columns_dos(raw, meta);

    % clean
    m = isfinite(E) & isfinite(DOS);
    if ~isempty(filling), m = m & isfinite(filling); end
    if ~isempty(doping),  m = m & isfinite(doping);  end

    E = E(m);
    DOS = DOS(m);
    if ~isempty(filling), filling = filling(m); end
    if ~isempty(doping),  doping  = doping(m);  end

    if isempty(DOS)
        error("No valid rows after cleaning: %s", in_path);
    end

    % -----------------------------
    % output dir
    % -----------------------------
    out_dir = "plot";
    if ~exist(out_dir,"dir"), mkdir(out_dir); end

    % -----------------------------
    % helper: title lines from header meta
    % -----------------------------
    function ttl = make_title(xname)
        lines = {};

        % line 1: model + xname
        if strlength(meta.model) > 0
            lines{end+1} = sprintf('%s DOS (%s)', upper(char(meta.model)), xname);
        else
            lines{end+1} = sprintf('DOS (%s)', xname);
        end

        % line 2: key parameters (if available)
        meta2 = {};
        if isfinite(meta.nLayer), meta2{end+1} = sprintf('Nl=%d', meta.nLayer); end %#ok<AGROW>
        if isfinite(meta.Dfield_eV), meta2{end+1} = sprintf('D=%.6g eV', meta.Dfield_eV); end %#ok<AGROW>
        if isfinite(meta.T_K), meta2{end+1} = sprintf('T=%.6g K', meta.T_K); end %#ok<AGROW>
        if strlength(meta.kmesh_type) > 0, meta2{end+1} = sprintf('kmesh=%s', char(meta.kmesh_type)); end %#ok<AGROW>
        if isfinite(meta.Nk), meta2{end+1} = sprintf('Nk=%d', meta.Nk); end %#ok<AGROW>
        if isfinite(meta.dk_frac), meta2{end+1} = sprintf('dk=%.6g', meta.dk_frac); end %#ok<AGROW>
        if isfinite(meta.eta), meta2{end+1} = sprintf('eta=%.6g', meta.eta); end %#ok<AGROW>
        if isfinite(meta.dE), meta2{end+1} = sprintf('dE=%.6g', meta.dE); end %#ok<AGROW>
        if ~isempty(meta2)
            lines{end+1} = strjoin(meta2, ', ');
        end

        % line 3: E_scan (if available)
        if isfinite(meta.E_low) && isfinite(meta.E_high) && isfinite(meta.num_e)
            lines{end+1} = sprintf('E\\_scan=[%.6g, %.6g], num\\_e=%d', meta.E_low, meta.E_high, meta.num_e);
        end

        ttl = strjoin(lines, newline);
    end

    FS = 14;

    % -----------------------------
    % 1) DOS vs E
    % -----------------------------
    [Es, idxE] = sort(E);
    DOSE = DOS(idxE);

    fig = figure('Color','w', 'Units','pixels', 'Position',[100 100 560 420]);
    ax = axes(fig); hold(ax,'on');

    plot(ax, Es, DOSE, 'LineWidth', 2.0);
    xlabel(ax, 'E (eV)', 'Interpreter','none', 'FontSize',FS);
    ylabel(ax, 'DOS (arb.)', 'Interpreter','none', 'FontSize',FS);
    title(ax, make_title('vs E'), 'Interpreter','none', 'FontSize',FS, 'FontWeight','normal');
    set(ax, 'FontSize',FS, 'LineWidth',1.0, 'TickDir','out', 'Box','on');
    grid(ax,'on'); xlim(ax,[min(Es) max(Es)]);

    out_png = fullfile(out_dir, base + "_vsE.png");
    exportgraphics(fig, out_png, 'Resolution',300);
    fprintf("Saved figure:\n  %s\n", out_png);
    out_png_list(end+1) = string(out_png);

    % -----------------------------
    % 2) DOS vs filling (if available)
    % -----------------------------
    if ~isempty(filling)
        [xs, idx] = sort(filling);
        ys = DOS(idx);

        fig = figure('Color','w', 'Units','pixels', 'Position',[120 120 560 420]);
        ax = axes(fig); hold(ax,'on');

        plot(ax, xs, ys, 'LineWidth', 2.0);
        xline(ax, 0.5, '--', 'LineWidth', 1.0); % neutrality convention
        xlabel(ax, 'filling', 'Interpreter','none', 'FontSize',FS);
        ylabel(ax, 'DOS (arb.)', 'Interpreter','none', 'FontSize',FS);
        title(ax, make_title('vs filling'), 'Interpreter','none', 'FontSize',FS, 'FontWeight','normal');
        set(ax, 'FontSize',FS, 'LineWidth',1.0, 'TickDir','out', 'Box','on');
        grid(ax,'on'); xlim(ax,[min(xs) max(xs)]);

        out_png = fullfile(out_dir, base + "_vsFilling.png");
        exportgraphics(fig, out_png, 'Resolution',300);
        fprintf("Saved figure:\n  %s\n", out_png);
        out_png_list(end+1) = string(out_png);
    else
        fprintf("NOTE: filling column not found -> skip DOS vs filling.\n");
    end

    % -----------------------------
    % 3) DOS vs doping (if available)
    % -----------------------------
    if ~isempty(doping)
        [xs, idx] = sort(doping);
        ys = DOS(idx);

        fig = figure('Color','w', 'Units','pixels', 'Position',[140 140 560 420]);
        ax = axes(fig); hold(ax,'on');

        plot(ax, xs, ys, 'LineWidth', 2.0);
        xline(ax, 0.0, '--', 'LineWidth', 1.0); % neutrality convention
        xlabel(ax, 'doping (1e12 cm^{-2})', 'Interpreter','none', 'FontSize',FS);
        ylabel(ax, 'DOS (arb.)', 'Interpreter','none', 'FontSize',FS);
        title(ax, make_title('vs doping'), 'Interpreter','none', 'FontSize',FS, 'FontWeight','normal');
        set(ax, 'FontSize',FS, 'LineWidth',1.0, 'TickDir','out', 'Box','on');
        grid(ax,'on'); xlim(ax,[min(xs) max(xs)]);

        out_png = fullfile(out_dir, base + "_vsDoping.png");
        exportgraphics(fig, out_png, 'Resolution',300);
        fprintf("Saved figure:\n  %s\n", out_png);
        out_png_list(end+1) = string(out_png);
    else
        fprintf("NOTE: doping column not found -> skip DOS vs doping.\n");
    end
end

% ============================================================
% parse header meta for DOS files
% - robust keyword-based parsing (NOT filename)
% ============================================================
function meta = parse_header_meta_dos(in_path)
    meta = struct();
    meta.model      = "";
    meta.nLayer     = NaN;
    meta.Dfield_eV  = NaN;
    meta.T_K        = NaN;
    meta.kmesh_type = "";
    meta.Nk         = NaN;
    meta.dk_frac    = NaN;

    meta.E_low  = NaN;
    meta.E_high = NaN;
    meta.num_e  = NaN;
    meta.eta    = NaN;

    meta.dE     = NaN;

    meta.columns_line = ""; % the full "# columns: ..." content (without '#')

    fid = fopen(in_path,'r');
    if fid < 0, error("Cannot open file: %s", in_path); end

    num = "([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)";

    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end
        s0 = strtrim(tline);
        if ~startsWith(s0, "#"), break; end

        s = strtrim(erase(s0, "#"));

        % model = sk
        if strlength(meta.model)==0 && contains(s, "model")
            tok = regexp(s, "(^|\s)model\s*=\s*([A-Za-z0-9_]+)", "tokens", "once");
            if ~isempty(tok), meta.model = string(tok{2}); end
        end

        % nLayer = 5  (from RG_Structure block)
        if ~isfinite(meta.nLayer) && contains(s, "nLayer")
            tok = regexp(s, "(^|\s)nLayer\s*=\s*(\d+)", "tokens", "once");
            if ~isempty(tok), meta.nLayer = str2double(tok{2}); end
        end

        % Dfield_eV = 0.067
        if ~isfinite(meta.Dfield_eV) && contains(s, "Dfield_eV")
            tok = regexp(s, "(^|\s)Dfield_eV\s*=\s*" + num, "tokens", "once");
            if ~isempty(tok), meta.Dfield_eV = str2double(tok{2}); end
        end

        % T_K = 0.4 K
        if ~isfinite(meta.T_K) && contains(s, "T_K")
            tok = regexp(s, "(^|\s)T_K\s*=\s*" + num, "tokens", "once");
            if ~isempty(tok), meta.T_K = str2double(tok{2}); end
        end

        % kmesh = BZ Nk=1000 dk_frac=0.004
        if (strlength(meta.kmesh_type)==0 || ~isfinite(meta.Nk) || ~isfinite(meta.dk_frac)) && contains(s, "kmesh")
            tok = regexp(s, "kmesh\s*=\s*([A-Za-z0-9_]+)", "tokens", "once");
            if ~isempty(tok), meta.kmesh_type = string(tok{1}); end

            tk = regexp(s, "Nk\s*=\s*(\d+)", "tokens", "once");
            if ~isempty(tk), meta.Nk = str2double(tk{1}); end

            td = regexp(s, "dk_frac\s*=\s*" + num, "tokens", "once");
            if ~isempty(td), meta.dk_frac = str2double(td{1}); end
        end

        % E_scan = [0.87, 0.92], num_e=501, eta=0.00015
        if ( ~isfinite(meta.E_low) || ~isfinite(meta.E_high) || ~isfinite(meta.num_e) || ~isfinite(meta.eta) ) ...
                && contains(s, "E_scan")
            tk = regexp(s, "E_scan\s*=\s*\[\s*" + num + "\s*,\s*" + num + "\s*\]", "tokens", "once");
            if ~isempty(tk)
                meta.E_low  = str2double(tk{1});
                meta.E_high = str2double(tk{2});
            end

            tn = regexp(s, "num_e\s*=\s*(\d+)", "tokens", "once");
            if ~isempty(tn), meta.num_e = str2double(tn{1}); end

            te = regexp(s, "eta\s*=\s*" + num, "tokens", "once");
            if ~isempty(te), meta.eta = str2double(te{1}); end
        end

        % dE = ...
        if ~isfinite(meta.dE) && startsWith(s, "dE")
            tdE = regexp(s, "dE\s*=\s*" + num, "tokens", "once");
            if ~isempty(tdE), meta.dE = str2double(tdE{1}); end
        end

        % columns: ...
        if strlength(meta.columns_line)==0 && startsWith(lower(s), "columns:")
            meta.columns_line = string(s);
        end
    end

    fclose(fid);
end

% ============================================================
% extract columns based on header columns line (preferred)
% fallback: infer by column count
% ============================================================
function [E, filling, doping, DOS] = extract_columns_dos(raw, meta)
    filling = [];
    doping  = [];

    ncol = size(raw,2);

    % Always: col1=i, col2=E
    E = raw(:,2);

    % ---------- header-driven ----------
    if strlength(meta.columns_line) > 0
        s = char(meta.columns_line);
        s = strrep(s, "columns:", "");
        s = strtrim(s);

        % normalize multiple spaces
        s = regexprep(s, "\s+", " ");

        % Tokenize by spaces
        toks = strsplit(s, " ");

        % Find indices (1-based within toks)
        % Example: i E(eV) filling doping_1e12cm^-2 DOS
        idx_fill = find(strcmpi(toks, "filling"), 1);
        idx_dop  = find(contains(lower(toks), "doping"), 1);
        idx_dos  = find(strcmpi(toks, "DOS") | strcmpi(toks, "dos"), 1);

        % Map token index -> column in raw:
        % toks{1} corresponds to raw(:,1), etc.
        if isempty(idx_dos)
            % if header exists but doesn't say DOS, fallback
            idx_dos = ncol;
        end

        if idx_dos <= ncol
            DOS = raw(:, idx_dos);
        else
            DOS = raw(:, ncol);
        end

        if ~isempty(idx_fill) && idx_fill <= ncol
            filling = raw(:, idx_fill);
        end
        if ~isempty(idx_dop) && idx_dop <= ncol
            doping = raw(:, idx_dop);
        end

        % sanity: if DOS is empty, fallback
        if isempty(DOS)
            DOS = raw(:, end);
        end
        return;
    end

    % ---------- fallback by ncol ----------
    if ncol == 3
        DOS = raw(:,3);
    elseif ncol == 4
        c3 = raw(:,3);
        c4 = raw(:,4);
        if median(abs(c3), 'omitnan') <= 2.0
            filling = c3;
            DOS = c4;
        else
            doping = c3;
            DOS = c4;
        end
    else
        filling = raw(:,3);
        doping  = raw(:,4);
        DOS     = raw(:,5);
    end
end
