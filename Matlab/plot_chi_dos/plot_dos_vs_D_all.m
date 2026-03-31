function out = plot_dos_vs_D_all()
% Plot DOS(E) curves for all D files in one folder, colored by D with colorbar.

    default_root = "/Users/haoranyan/rg_master/data/";
    if ~isfolder(default_root)
        default_root = pwd;
    end

    root = uigetdir(default_root, 'Select folder containing DOS files');
    if isequal(root,0)
        error('No folder selected.');
    end

    files = dir(fullfile(root, 'dos_D*meV*.txt'));
    if isempty(files)
        error('No files matching dos_D*meV*.txt found in:\n%s', root);
    end

    % ------------------------------------------------------------
    % collect all curves
    % ------------------------------------------------------------
    S = struct('name', {}, 'path', {}, 'D_meV', {}, 'E', {}, 'DOS', {});
    nkeep = 0;

    for i = 1:numel(files)
        fname = files(i).name;
        fpath = fullfile(files(i).folder, fname);

        D_meV = parse_D_from_filename_(fname);
        if isnan(D_meV)
            fprintf('[skip] cannot parse D from file name: %s\n', fname);
            continue;
        end

        M = read_numeric_block_(fpath);
        if isempty(M) || size(M,2) < 2
            fprintf('[skip] no valid numeric block: %s\n', fname);
            continue;
        end

        E   = M(:,2);
        DOS = M(:,end);

        good = isfinite(E) & isfinite(DOS);
        E = E(good);
        DOS = DOS(good);

        if isempty(E)
            fprintf('[skip] empty valid E/DOS: %s\n', fname);
            continue;
        end

        nkeep = nkeep + 1;
        S(nkeep).name  = fname;
        S(nkeep).path  = fpath;
        S(nkeep).D_meV = D_meV;
        S(nkeep).E     = E;
        S(nkeep).DOS   = DOS;
    end

    if isempty(S)
        error('No valid DOS files could be read.');
    end

    % ------------------------------------------------------------
    % sort by D
    % ------------------------------------------------------------
    [~, ord] = sort([S.D_meV]);
    S = S(ord);
    Dvals = [S.D_meV];

    % ------------------------------------------------------------
    % plot all in one figure with colorbar
    % ------------------------------------------------------------
    fig = figure('Color', 'w');
    ax = axes(fig);
    hold(ax, 'on');

    cmap = jet(256);
    colormap(ax, cmap);

    Dmin = min(Dvals);
    Dmax = max(Dvals);

    if abs(Dmax - Dmin) < 1e-12
        Dnorm = 0.5 * ones(size(Dvals));
    else
        Dnorm = (Dvals - Dmin) / (Dmax - Dmin);
    end

    for i = 1:numel(S)
        idx = max(1, min(size(cmap,1), 1 + round(Dnorm(i) * (size(cmap,1)-1))));
        thisColor = cmap(idx, :);

        plot(ax, S(i).E, S(i).DOS, ...
            'LineWidth', 1.5, ...
            'Color', thisColor);
    end

    xlabel(ax, 'E (eV)', 'FontSize', 16);
    ylabel(ax, 'DOS', 'FontSize', 16);
    title(ax, 'DOS for different D', 'FontSize', 16);
    set(ax, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
    grid(ax, 'on');

    caxis(ax, [Dmin, Dmax]);
    cb = colorbar(ax);
    cb.Label.String = 'D (meV)';
    cb.FontSize = 14;


end


function D_meV = parse_D_from_filename_(fname)
    tok = regexp(fname, 'dos_D([+-]?\d+(?:\.\d+)?)meV', 'tokens', 'once');
    if isempty(tok)
        D_meV = NaN;
    else
        D_meV = str2double(tok{1});
    end
end


function M = read_numeric_block_(fpath)
    fid = fopen(fpath, 'r');
    if fid < 0
        error('Cannot open file: %s', fpath);
    end

    C = {};
    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end

        s = strtrim(tline);
        if isempty(s), continue; end
        if startsWith(s, '#'), continue; end

        nums = sscanf(s, '%f').';
        if ~isempty(nums)
            C{end+1,1} = nums; %#ok<AGROW>
        end
    end
    fclose(fid);

    if isempty(C)
        M = [];
        return;
    end

    ncol = max(cellfun(@numel, C));
    keep = cellfun(@numel, C) == ncol;
    C = C(keep);

    if isempty(C)
        M = [];
        return;
    end

    M = vertcat(C{:});
end