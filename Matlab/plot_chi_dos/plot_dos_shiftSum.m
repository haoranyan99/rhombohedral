function out_png_list = plot_dos_shiftMean_multiShift_vs_x()
% Plot multiple shift-mean DOS curves vs x-axis = doping or mu
%
% Rule:
%   1) shift is always done on the ORIGINAL energy axis E
%   2) for each target E_i, evaluate left/right points at E_i +/- dE
%   3) after shifting in energy, parse x = doping or mu at shifted points
%   4) if exact shifted energy is not found, use mean of nearest neighbors
%
% Output:
%   x_plot(i)   = 0.5 * [ x(E_i + dE) + x(E_i - dE) ]
%   y_plot(i)   = 0.5 * [ DOS(E_i + dE) + DOS(E_i - dE) ]
%
% New features:
%   - support multiple shift values
%   - use hot colormap
%   - shift = 0 corresponds to the WHITE end of the colormap
%   - but the actual 0-shift curve is plotted as gray dashed line
%     representing the original DOS

out_png_list = strings(0);

% ==========================
% USER SETTINGS
% ==========================
start_dir = "E:/rg_master/data/";
mu_shift_list_meV = [0, 2, 4, 6, 8, 10];   % multiple shifts in meV
save_png = false;
FS = 14;
LW_org = 2.0;
LW_shift = 2.2;

% ==========================
% Choose x-axis by UI
% ==========================
x_mode = pick_xmode_();   % "doping" or "mu"
fprintf("[x_mode] %s\n", x_mode);

% sort + unique
mu_shift_list_meV = unique(mu_shift_list_meV(:).', 'sorted');

fprintf("[shift list meV] ");
fprintf("%.6g ", mu_shift_list_meV);
fprintf("\n");

% ==========================
% Select files
% ==========================
[file, path] = uigetfile("*.txt", "Select DOS txt file(s)", start_dir, "MultiSelect","on");
if isequal(file,0)
    error("No file selected.");
end
if ischar(file)
    file = {file};
end

% ==========================
% Prepare colormap
% shift = 0 should be white end
% hot() default is black -> red -> yellow -> white
% so we flip it to get white -> yellow -> red -> black
% ==========================
shift_min = min(mu_shift_list_meV);
shift_max = max(mu_shift_list_meV);

n_cmap = 256;
cmap_full = flipud(hot(n_cmap));

if shift_min == shift_max
    get_shift_color = @(s) cmap_full(1,:);
else
    get_shift_color = @(s) interp1( ...
        linspace(shift_min, shift_max, n_cmap), ...
        cmap_full, ...
        s, ...
        'linear', 'extrap');
end

% ==========================
% Loop over files
% ==========================
for k = 1:numel(file)
    in_path = fullfile(path, file{k});
    [~, base, ~] = fileparts(in_path);

    raw = readmatrix(in_path, "FileType","text", "CommentStyle","#");
    if isempty(raw) || size(raw,2) < 3
        fprintf("[skip] bad file: %s\n", in_path);
        continue;
    end

    cols_line = parse_columns_line_(in_path);
    meta_dos  = parse_header_mu_doping_(in_path);

    try
        [E, DOS, x_raw] = extract_E_DOS_x_(raw, cols_line, x_mode, meta_dos);
    catch ME
        fprintf("[skip] %s | %s\n", in_path, ME.message);
        continue;
    end

    m = isfinite(E) & isfinite(DOS) & isfinite(x_raw);
    E     = E(m);
    DOS   = DOS(m);
    x_raw = x_raw(m);

    if numel(E) < 3
        fprintf("[skip] too few valid points: %s\n", in_path);
        continue;
    end

    % sort by energy
    [E, ord] = sort(E);
    DOS   = DOS(ord);
    x_raw = x_raw(ord);

    % collapse repeated E by averaging DOS and x on same E
    [E, DOS, x_raw] = average_by_E_(E, DOS, x_raw);

    if numel(E) < 3
        fprintf("[skip] too few unique E points: %s\n", in_path);
        continue;
    end

    % original curve
    [x_org, y_org] = average_by_x_(x_raw, DOS);

    % ==========================
    % Plot
    % ==========================
    fig = figure('Color','w','Units','pixels','Position',[120 120 860 580], ...
        'Name', sprintf('multi-shift DOS vs %s', x_mode));

    ax = axes(fig);
    hold(ax,'on');
    box(ax,'on');
    grid(ax,'on');
    set(ax,'FontSize',FS,'LineWidth',1,'TickDir','out');

    % plot original DOS if shift=0 is included
    if any(mu_shift_list_meV == 0)
        plot(ax, x_org, y_org, '--', ...
            'Color', [0.45 0.45 0.45], ...
            'LineWidth', LW_org, ...
            'DisplayName', 'original DOS (\DeltaE = 0)');
    end

    % plot nonzero shifted curves
    for is = 1:numel(mu_shift_list_meV)
        shift_meV = mu_shift_list_meV(is);

        if shift_meV == 0
            continue;
        end

        this_color = get_shift_color(shift_meV);
        shift_eV = 1e-3 * shift_meV;

        x_plot = nan(size(E));
        y_plot = nan(size(E));

        for i = 1:numel(E)
            E0 = E(i);

            [DOSp, xp, okp] = sample_at_shifted_energy_(E, DOS, x_raw, E0 + shift_eV);
            [DOSm, xm, okm] = sample_at_shifted_energy_(E, DOS, x_raw, E0 - shift_eV);

            if okp && okm
                y_plot(i) = 0.5 * (DOSp + DOSm);
                x_plot(i) = 0.5 * (xp   + xm);
            end
        end

        good = isfinite(x_plot) & isfinite(y_plot);
        x_plot = x_plot(good);
        y_plot = y_plot(good);

        if isempty(x_plot)
            fprintf("[skip] no valid shifted points: %s | shift=%.6g meV\n", in_path, shift_meV);
            continue;
        end

        % average repeated x after shift
        [x_plot, y_plot] = average_by_x_(x_plot, y_plot);

        plot(ax, x_plot, y_plot, '-', ...
            'Color', this_color, ...
            'LineWidth', LW_shift, ...
            'DisplayName', sprintf('\\DeltaE = %.4g meV', shift_meV));
    end

    % labels
    if x_mode == "doping"
        xlabel(ax, 'doping (1e12 cm^{-2})', 'Interpreter','tex', 'FontSize',FS);
    else
        xlabel(ax, '\mu (eV)', 'Interpreter','tex', 'FontSize',FS);
    end

    ylabel(ax, 'DOS', 'FontSize',FS);
    title(ax, sprintf('%s | multi-shift mean DOS vs %s', base, x_mode), ...
        'Interpreter','none', 'FontWeight','normal');

    % colorbar
    colormap(ax, cmap_full);
    cb = colorbar(ax);
    cb.Label.String = '\DeltaE (meV)';
    cb.FontSize = FS;

    if shift_min == shift_max
        caxis(ax, [shift_min - 1, shift_max + 1]);
    else
        caxis(ax, [shift_min, shift_max]);
    end

    cb.Ticks = mu_shift_list_meV;

    legend(ax, 'Location', 'best');

    if save_png
        out_png = fullfile(path, base + "_multiShiftMean_vs_" + x_mode + ".png");
        exportgraphics(fig, out_png, 'Resolution', 300);
        out_png_list(end+1) = string(out_png);
    end
end

end

% ============================================================
% helpers
% ============================================================

function x_mode = pick_xmode_()
c = questdlg('Choose x-axis:', 'x-axis', 'doping','mu','doping');
if isempty(c)
    x_mode = "doping";
else
    x_mode = string(c);
end
end

function cols_line = parse_columns_line_(in_path)
cols_line = "";
fid = fopen(in_path,'r');
if fid < 0
    return;
end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    s = strtrim(ln);
    if ~startsWith(s,"#"), break; end
    s2 = strtrim(erase(s,"#"));
    if startsWith(lower(s2), "columns:")
        cols_line = string(s2);
        return;
    end
end
end

function meta = parse_header_mu_doping_(fpath)
meta = struct('doping',NaN,'mu',NaN);

fid = fopen(fpath,'r');
if fid < 0
    return;
end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

num = '([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)';

for t = 1:4000
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    s = strtrim(ln);
    if isempty(s), continue; end
    if s(1) ~= '#', break; end

    if ~isfinite(meta.doping)
        tok = regexp(s, ['^\#\s*doping\s*=\s*' num], 'tokens','once');
        if ~isempty(tok)
            meta.doping = str2double(tok{1});
        end
    end

    if ~isfinite(meta.mu)
        tok = regexp(s, ['^\#\s*mu\s*=\s*' num], 'tokens','once');
        if isempty(tok)
            tok = regexp(s, ['^\#\s*EF\s*=\s*' num], 'tokens','once');
        end
        if ~isempty(tok)
            meta.mu = str2double(tok{1});
        end
    end

    if isfinite(meta.doping) && isfinite(meta.mu)
        break;
    end
end
end

function [E, DOS, x] = extract_E_DOS_x_(raw, cols_line, x_mode, meta_dos)
% Parse:
%   E   from numeric column if possible
%   DOS from DOS column if possible
%   x   = doping or mu, numeric column preferred, header fallback
%
% Common numeric format:
%   i  E  filling/mu  doping  DOS

ncol = size(raw,2);

E   = raw(:,2);
DOS = raw(:,end);
x   = nan(size(E));

if strlength(cols_line) > 0
    s = char(cols_line);
    s = strrep(s, "columns:", "");
    s = strtrim(s);
    s = regexprep(s, "\s+", " ");
    toks = strsplit(s, " ");

    idx_E = find(strcmpi(toks,"E") | contains(lower(toks),"e(ev)"), 1);
    if ~isempty(idx_E) && idx_E <= ncol
        E = raw(:, idx_E);
    end

    idx_DOS = find(strcmpi(toks,"DOS") | strcmpi(toks,"dos"), 1);
    if isempty(idx_DOS)
        idx_DOS = ncol;
    end
    DOS = raw(:, min(idx_DOS, ncol));

    if x_mode == "doping"
        idx_x = find(contains(lower(toks),"doping"), 1);
        if ~isempty(idx_x) && idx_x <= ncol
            x = raw(:, idx_x);
            return;
        end
        if isfinite(meta_dos.doping)
            x = meta_dos.doping * ones(size(E));
            return;
        end
        error("cannot find doping column in numeric table nor 'doping=' in header.");
    else
        idx_x = find(strcmpi(toks,"mu") | strcmpi(toks,"EF") | ...
                     contains(lower(toks),"mu") | contains(lower(toks),"ef") | ...
                     contains(lower(toks),"filling"), 1);
        if ~isempty(idx_x) && idx_x <= ncol
            x = raw(:, idx_x);
            return;
        end
        if isfinite(meta_dos.mu)
            x = meta_dos.mu * ones(size(E));
            return;
        end
        error("cannot find mu/EF/filling column in numeric table nor 'mu='/'EF=' in header.");
    end
end

% no columns line: fallback
if x_mode == "doping"
    if ncol >= 5
        E   = raw(:,2);
        x   = raw(:,4);
        DOS = raw(:,5);
        return;
    elseif ncol == 4
        E   = raw(:,2);
        x   = raw(:,3);
        DOS = raw(:,4);
        return;
    elseif isfinite(meta_dos.doping)
        E   = raw(:,2);
        DOS = raw(:,end);
        x   = meta_dos.doping * ones(size(E));
        return;
    else
        error("no '# columns:' and cannot infer doping column; also no header doping.");
    end
else
    if ncol >= 5
        E   = raw(:,2);
        x   = raw(:,3);
        DOS = raw(:,5);
        return;
    elseif isfinite(meta_dos.mu)
        E   = raw(:,2);
        DOS = raw(:,end);
        x   = meta_dos.mu * ones(size(E));
        return;
    else
        error("no '# columns:' and cannot infer mu column; also no header mu/EF.");
    end
end
end

function [E_u, DOS_u, x_u] = average_by_E_(E, DOS, x)
Eq = round(E, 12);
[Es, ord] = sort(Eq);
DOSs = DOS(ord);
xs   = x(ord);

[Eu, ~, ic] = unique(Es);
DOS_u = accumarray(ic, DOSs, [], @(v) mean(v,'omitnan'));
x_u   = accumarray(ic, xs,   [], @(v) mean(v,'omitnan'));
E_u   = Eu;
end

function [x_u, y_u] = average_by_x_(x, y)
if isempty(x)
    x_u = x;
    y_u = y;
    return;
end

xq = round(x, 12);
[xs, ord] = sort(xq);
ys = y(ord);

[ux, ~, ic] = unique(xs);
y_u = accumarray(ic, ys, [], @(v) mean(v,'omitnan'));
x_u = ux;
end

function [DOSv, xv, ok] = sample_at_shifted_energy_(E, DOS, x, Etarget)
% Exact hit preferred.
% If not found, use mean of nearest two neighbors around Etarget.
% If only one side exists, use nearest single point.

ok = false;
DOSv = NaN;
xv   = NaN;

if isempty(E)
    return;
end

% exact / very close hit
tol = max(1e-12, 1e-9 * max(1, max(abs(E))));
id_exact = find(abs(E - Etarget) <= tol, 1, 'first');
if ~isempty(id_exact)
    DOSv = DOS(id_exact);
    xv   = x(id_exact);
    ok = isfinite(DOSv) && isfinite(xv);
    return;
end

% neighbors
iR = find(E > Etarget, 1, 'first');
iL = find(E < Etarget, 1, 'last');

if ~isempty(iL) && ~isempty(iR)
    DOSv = mean([DOS(iL), DOS(iR)], 'omitnan');
    xv   = mean([x(iL),   x(iR)],   'omitnan');
    ok = isfinite(DOSv) && isfinite(xv);
    return;
elseif ~isempty(iL)
    DOSv = DOS(iL);
    xv   = x(iL);
    ok = isfinite(DOSv) && isfinite(xv);
    return;
elseif ~isempty(iR)
    DOSv = DOS(iR);
    xv   = x(iR);
    ok = isfinite(DOSv) && isfinite(xv);
    return;
end
end