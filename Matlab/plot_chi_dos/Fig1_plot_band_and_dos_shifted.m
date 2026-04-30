function Fig1_plot_band_and_dos_shifted()

clc; close all;

% =============================
% User settings
% =============================
default_root = "E:/rg_master/data/";
if ~isfolder(default_root)
    default_root = "/Users/haoranyan/rg_master/data/";
end
if ~isfolder(default_root)
    default_root = pwd;
end

FS = 14;
cut_K_to_Kp_half = true;

% =============================
% Select DOS file
% =============================
[fname_dos, fpath_dos] = uigetfile({'*.txt','DOS data (*.txt)'}, ...
    'Select DOS file', fullfile(default_root, '*.txt'));

if isequal(fname_dos,0)
    return;
end

dos_path = fullfile(fpath_dos, fname_dos);

% =============================
% Select band file
% =============================
[fname_band, fpath_band] = uigetfile({'*.txt','Band data (*.txt)'}, ...
    'Select bare band file used for mu shift', fullfile(default_root, '*.txt'));

if isequal(fname_band,0)
    return;
end

band_path = fullfile(fpath_band, fname_band);

% =============================
% Read band and compute mu shift
% =============================
A_band = readmatrix(band_path, 'CommentStyle', '#');

if isempty(A_band) || size(A_band,2) < 11
    error('Band file must contain columns: ik s kx ky E0 E1 ... at least E6.');
end

s_band = A_band(:,2);
E_band = A_band(:,5:end);
norb = size(E_band,2);

if norb < 7
    error('Need at least 7 bands: E0...E6.');
end

% File bands are E0,E1,... so E5 -> column 6, E6 -> column 7
E5_max = max(E_band(:,6), [], 'omitnan');
E6_min = min(E_band(:,7), [], 'omitnan');

mu_shift = 0.5 * (E5_max + E6_min);
E_band = E_band - mu_shift;

fprintf('mu_shift = %.10f eV\n', mu_shift);
fprintf('After shift: max(E5)=%.10f eV, min(E6)=%.10f eV\n', ...
    E5_max - mu_shift, E6_min - mu_shift);

% =============================
% Parse band xticks
% =============================
[xtick_pos, xtick_lab] = parse_xticks_band(band_path, s_band);

xtick_pos(1) = s_band(1);
xtick_pos(end) = s_band(end);

% =============================
% Cut K -> K' half
% =============================
if cut_K_to_Kp_half && numel(xtick_pos) >= 3
    K_pos  = xtick_pos(2);
    Kp_pos = xtick_pos(3);
    cut_pos = 0.5 * (K_pos + Kp_pos);

    mask_band = s_band <= cut_pos;
    s_band = s_band(mask_band);
    E_band = E_band(mask_band,:);

    xtick_pos = [xtick_pos(1), K_pos, cut_pos];
    xtick_lab = {'M','K',"mid(K,K')"};
end

% =============================
% Read DOS
% =============================
meta = parse_header_meta_dos(dos_path);

raw = readmatrix(dos_path, "FileType","text", "CommentStyle","#");
if isempty(raw) || size(raw,2) < 3
    error("Bad DOS format: need at least 3 numeric columns.");
end

[E_dos, filling, doping, DOS] = extract_columns_dos(raw, meta);

mask = isfinite(E_dos) & isfinite(DOS);
if ~isempty(filling)
    mask = mask & isfinite(filling);
end
if ~isempty(doping)
    mask = mask & isfinite(doping);
end

E_dos = E_dos(mask);
DOS = DOS(mask);

E_dos = E_dos - mu_shift;

% =============================
% Sort DOS
% =============================
[Es, idx] = sort(E_dos);
DOSE = DOS(idx);

dos_xlim = [min(Es), max(Es)];

% =============================
% Plot 1: shifted DOS
% =============================
fig1 = figure('Color','w','Units','pixels','Position',[100 100 620 460]);
ax1 = axes(fig1);
hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');

plot(ax1, Es, DOSE, 'LineWidth', 2.0);
xline(ax1, 0, '--', 'LineWidth', 1.2);

xlabel(ax1, 'E - \mu (eV)', 'FontSize', FS);
ylabel(ax1, 'DOS (arb.)', 'FontSize', FS);
title(ax1, sprintf('Shifted DOS, \\mu = %.6f eV', mu_shift), ...
    'FontSize', FS, 'FontWeight','normal');

set(ax1, ...
    'FontSize', FS, ...
    'LineWidth', 1.1, ...
    'TickDir','out', ...
    'Box','on');

xlim(ax1, dos_xlim);

% =============================
% Plot 2: shifted band
% =============================
fig2 = figure('Color','w','Units','pixels','Position',[160 120 820 520]);
ax2 = axes(fig2);
hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');

for ib = 1:norb
    plot(ax2, s_band, E_band(:,ib), 'LineWidth', 1.2);
end

yline(ax2, 0, '--', 'LineWidth', 1.2);

xlabel(ax2, 'k-line coordinate s (1/\AA)', 'FontSize', FS);
ylabel(ax2, 'E - \mu (eV)', 'FontSize', FS);
title(ax2, sprintf('Shifted bare band, \\mu = %.6f eV', mu_shift), ...
    'FontSize', FS, 'FontWeight','normal');

set(ax2, ...
    'FontSize', FS, ...
    'LineWidth', 1.1, ...
    'TickDir','out', ...
    'Box','on', ...
    'XTick', xtick_pos, ...
    'XTickLabel', xtick_lab, ...
    'XLim', [s_band(1), s_band(end)]);

ylim(ax2, dos_xlim);

for i = 2:numel(xtick_pos)-1
    xline(ax2, xtick_pos(i), '-', 'Color', [0.75 0.75 0.75], ...
        'LineWidth', 0.8, 'HandleVisibility','off');
end

end

% ============================================================
% Helper: parse band xticks
% ============================================================
function [xtick_pos, xtick_lab] = parse_xticks_band(in_path, s)

xtick_pos = [];
xtick_lab = {};

fid = fopen(in_path, 'r');
if fid < 0
    error('Cannot open file: %s', in_path);
end

while ~feof(fid)
    line = strtrim(fgetl(fid));

    if startsWith(line, '# xticks:')
        blocks = regexp(line, '\(([^)]+)\)', 'tokens');

        for i = 1:numel(blocks)
            item = strtrim(blocks{i}{1});
            parts = strsplit(item, ',');

            if numel(parts) >= 2
                lab = strtrim(parts{1});
                pos = str2double(strtrim(parts{2}));

                if ~isnan(pos)
                    xtick_lab{end+1} = lab;
                    xtick_pos(end+1) = pos;
                end
            end
        end
        break;
    end
end

fclose(fid);

if isempty(xtick_pos)
    xtick_pos = linspace(s(1), s(end), 3);
    xtick_lab = {'start','mid','end'};
end

end

% ============================================================
% Helper: parse DOS header
% ============================================================
function meta = parse_header_meta_dos(in_path)

meta = struct();
meta.columns_line = "";

fid = fopen(in_path,'r');
if fid < 0
    error("Cannot open file: %s", in_path);
end

while true
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end

    s0 = strtrim(tline);
    if ~startsWith(s0, "#")
        break;
    end

    s = strtrim(erase(s0, "#"));

    if strlength(meta.columns_line)==0 && startsWith(lower(s), "columns:")
        meta.columns_line = string(s);
    end
end

fclose(fid);

end

% ============================================================
% Helper: extract DOS columns
% ============================================================
function [E, filling, doping, DOS] = extract_columns_dos(raw, meta)

filling = [];
doping  = [];
ncol = size(raw,2);

E = raw(:,2);

if isfield(meta, "columns_line") && strlength(meta.columns_line) > 0
    s = char(meta.columns_line);
    s = strrep(s, "columns:", "");
    s = strtrim(s);
    s = regexprep(s, "\s+", " ");
    toks = strsplit(s, " ");

    idx_fill = find(strcmpi(toks, "filling"), 1);
    idx_dop  = find(contains(lower(toks), "doping"), 1);
    idx_dos  = find(strcmpi(toks, "DOS") | strcmpi(toks, "dos"), 1);

    if isempty(idx_dos)
        idx_dos = ncol;
    end

    DOS = raw(:, min(idx_dos, ncol));

    if ~isempty(idx_fill) && idx_fill <= ncol
        filling = raw(:, idx_fill);
    end

    if ~isempty(idx_dop) && idx_dop <= ncol
        doping = raw(:, idx_dop);
    end

    return;
end

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