function out = plot_dos_vs_D_all()
% Plot DOS(E) curves for selected D range, and only plot N curves.
% All parameters are set inside this script, no input arguments needed.

    % ------------------------------------------------------------
    % default settings
    % ------------------------------------------------------------
    default_root = "/Users/haoranyan/rg_master/data/";
    D_range = [-300, 0];   % only keep D in this range (meV)
    N_plot  = 8;           % only plot N curves within the range

    % if true, ask folder by uigetdir; otherwise use default_root directly
    use_uigetdir = true;

    % ------------------------------------------------------------
    % choose root folder
    % ------------------------------------------------------------
    if ~isfolder(default_root)
        default_root = pwd;
    end

    if use_uigetdir
        root = uigetdir(default_root, 'Select folder containing DOS files');
        if isequal(root,0)
            error('No folder selected.');
        end
    else
        root = default_root;
    end

    files = dir(fullfile(root, 'dos_D*meV*.txt'));
    if isempty(files)
        error('No files matching dos_D*meV*.txt found in:\n%s', root);
    end

    % ------------------------------------------------------------
    % collect all curves in the D range
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

        % keep only D in range
        if D_meV < D_range(1) || D_meV > D_range(2)
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
        error('No valid DOS files found in D range [%g, %g] meV.', ...
            D_range(1), D_range(2));
    end

    % ------------------------------------------------------------
    % sort by D
    % ------------------------------------------------------------
    [~, ord] = sort([S.D_meV]);
    S = S(ord);
    Dvals_all = [S.D_meV];

    % ------------------------------------------------------------
    % choose only N curves uniformly from the range
    % ------------------------------------------------------------
    nall = numel(S);

    if N_plot >= nall
        pick_idx = 1:nall;
    else
        pick_idx = round(linspace(1, nall, N_plot));
        pick_idx = unique(pick_idx, 'stable');

        % if rounding caused fewer than N points, fill in
        if numel(pick_idx) < N_plot
            rest = setdiff(1:nall, pick_idx, 'stable');
            nneed = N_plot - numel(pick_idx);
            pick_idx = [pick_idx, rest(1:nneed)];
            pick_idx = sort(pick_idx);
        end
    end

    Splot = S(pick_idx);
    Dvals = [Splot.D_meV];

    fprintf('Total valid curves in range [%g, %g] meV: %d\n', ...
        D_range(1), D_range(2), nall);
    fprintf('Number of curves plotted: %d\n', numel(Splot));
    fprintf('Plotted D values (meV):\n');
    disp(Dvals);

    % ------------------------------------------------------------
    % plot
    % ------------------------------------------------------------
    fig = figure('Color', 'w');
    ax = axes(fig);
    hold(ax, 'on');

    cmap = hot(256);
    colormap(ax, cmap);

    Dmin = min(Dvals);
    Dmax = max(Dvals);

    if abs(Dmax - Dmin) < 1e-12
        Dnorm = 0.5 * ones(size(Dvals));
    else
        Dnorm = (Dvals - Dmin) / (Dmax - Dmin);
    end

    for i = 1:numel(Splot)
        idx = max(1, min(size(cmap,1), 1 + round(Dnorm(i) * (size(cmap,1)-1))));
        thisColor = cmap(idx, :);

        plot(ax, Splot(i).E, Splot(i).DOS, ...
            'LineWidth', 1.5, ...
            'Color', thisColor, ...
            'DisplayName', sprintf('D = %.3f meV', Splot(i).D_meV));
    end

    xlabel(ax, 'E (eV)', 'FontSize', 16);
    ylabel(ax, 'DOS', 'FontSize', 16);
    title(ax, sprintf('DOS for D in [%g, %g] meV (%d curves)', ...
        D_range(1), D_range(2), numel(Splot)), 'FontSize', 16);

    set(ax, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
    grid(ax, 'on');

    caxis(ax, [Dmin, Dmax]);
    cb = colorbar(ax);
    cb.Label.String = 'D (meV)';
    cb.FontSize = 14;

    % optional legend
    % legend(ax, 'show', 'Location', 'best');

    % ------------------------------------------------------------
    % output
    % ------------------------------------------------------------
    out.root      = root;
    out.D_range   = D_range;
    out.N_plot    = N_plot;
    out.S_all     = S;
    out.S_plot    = Splot;
    out.pick_idx  = pick_idx;
    out.Dvals_all = Dvals_all;
    out.Dvals     = Dvals;
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