function plot_test_rg_structure(nLayer, pressure)
% plot_test_rg_structure
% Read single dump file:
%   data/test_Structure_and_KPath_dump.txt
% Sections used:
%   STRUCTURE_META, ATOMS, REAL_CELL
%
% Inputs:
%   nLayer   : number of layers (for filename only; if omitted, infer from data/meta)
%   pressure : GPa (for filename only; if omitted, infer from meta; otherwise 0)

    if nargin < 1 || isempty(nLayer)
        nLayer = [];
    end
    if nargin < 2 || isempty(pressure)
        pressure = [];
    end

    % ===============================
    % Font size settings (EDIT HERE)
    % ===============================
    fs_axis   = 16;   % tick labels
    fs_label  = 16;   % xlabel / ylabel
    fs_title  = 16;   % title
    fs_atom   = 12;   % atom text "C"

    dumpfile = "data/test_Structure_and_KPath_dump.txt";
    if ~isfile(dumpfile)
        error("Cannot find dump file: %s", dumpfile);
    end

    % -------- parse sections --------
    meta = read_section_numeric(dumpfile, "STRUCTURE_META"); % 1 row: nLayer pressure a d vacuum
    atoms = read_section_numeric(dumpfile, "ATOMS");         % columns: x y z layer sub
    cellv = read_section_numeric(dumpfile, "REAL_CELL");     % 1 row: a1x a1y a2x a2y

    if isempty(atoms) || size(atoms,2) < 5
        error("ATOMS section is empty or has wrong columns.");
    end
    if isempty(cellv) || numel(cellv) < 4
        error("REAL_CELL section is empty or has wrong columns.");
    end

    x = atoms(:,1);
    y = atoms(:,2);
    layer = atoms(:,4);

    % infer nLayer/pressure if not provided
    if isempty(nLayer)
        if ~isempty(meta) && numel(meta) >= 1
            nLayer = round(meta(1));
        else
            nLayer = max(layer) + 1;
        end
    end
    if isempty(pressure)
        if ~isempty(meta) && numel(meta) >= 2
            pressure = meta(2);
        else
            pressure = 0.0;
        end
    end

    a1 = cellv(1,1:2);
    a2 = cellv(1,3:4);

    % ---- ensure plot folder ----
    if ~exist("plot", "dir")
        mkdir("plot");
    end

    % ---- figure layout: 1 x num_layers ----
    fig = figure("Color","w");
    tiledlayout(1, nLayer, ...
        "Padding","compact", ...
        "TileSpacing","compact");

    % sharex/sharey: global limits
    pad = 0.8;
    xlim_all = [min(x)-pad, max(x)+pad];
    ylim_all = [min(y)-pad, max(y)+pad];

    for i = 0:nLayer-1
        ax = nexttile;
        hold(ax,"on");
        box(ax,"on");
        axis(ax,"equal");
        set(ax, "FontSize", fs_axis);

        % unit cell polygon
        origin = [0,0];
        cell_points = [ origin;
                        origin + a1;
                        origin + a1 + a2;
                        origin + a2;
                        origin ];
        plot(ax, cell_points(:,1), cell_points(:,2), "--k", "LineWidth", 1.5);

        % atoms in this layer
        idx = (layer == i);
        xi = x(idx);
        yi = y(idx);

        scatter(ax, xi, yi, 80, "filled");

        % atom labels
        for j = 1:numel(xi)
            text(ax, xi(j) + 0.1, yi(j) + 0.1, "C", "FontSize", fs_atom);
        end

        title(ax, sprintf("Layer %d", i+1), ...
            "FontSize", fs_title, "FontWeight","normal");

        xlabel(ax, "x (\AA)", "FontSize", fs_label);
        if i == 0
            ylabel(ax, "y (\AA)", "FontSize", fs_label);
        end

        xlim(ax, xlim_all);
        ylim(ax, ylim_all);
    end

    % ---- save ----
    outname = sprintf("./plot/test_RG_%dL_layers_%.3fGPa.png", nLayer, pressure);
    exportgraphics(fig, outname, "Resolution", 300);
    fprintf("Saved figure: %s\n", outname);
end

% ============================================================
% Local helper: read one numeric section from the dump file
% ============================================================
function A = read_section_numeric(filename, section_name)
% Reads numeric rows from:
%   # === SECTION_NAME ===
% until next "# ===" or EOF.
% Ignores comment lines starting with '#'.
% Returns numeric matrix A (can be empty).

    lines = readlines(filename);

    key = "# === " + string(section_name) + " ===";
    idx0 = find(lines == key, 1, "first");
    if isempty(idx0)
        A = [];
        return;
    end

    % find end marker (next section)
    idx1 = find(startsWith(lines(idx0+1:end), "# === "), 1, "first");
    if isempty(idx1)
        block = lines(idx0+1:end);
    else
        block = lines(idx0+1 : idx0+idx1-1);
    end

    % remove empty and comment lines
    block = strip(block);
    block = block(block ~= "");
    block = block(~startsWith(block, "#"));
    if isempty(block)
        A = [];
        return;
    end

    % parse numeric rows
    % each line -> sscanf -> row vector
    rows = cell(numel(block), 1);
    ncol = 0;
    for i = 1:numel(block)
        v = sscanf(block(i), "%f").';
        rows{i} = v;
        ncol = max(ncol, numel(v));
    end

    % build matrix (pad with NaN if inconsistent)
    A = nan(numel(rows), ncol);
    for i = 1:numel(rows)
        v = rows{i};
        A(i,1:numel(v)) = v;
    end

    % If only one row, keep as 1xN (fine)
end
