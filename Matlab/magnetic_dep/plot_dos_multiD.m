function out = plot_dos_multiD()
% Plot DOS vs doping or mu for selected D range, and only plot N curves.
% Before plotting, a dialog asks whether to use doping or mu as x-axis.

    % ------------------------------------------------------------
    % default settings
    % ------------------------------------------------------------
    default_root = "/Users/haoranyan/rg_master/data/";
    D_range = [-300, 0];   % only keep D in this range (meV)
    N_plot  = 8;           % only plot N curves within the range

    % optional x-range restriction
    use_x_range = false;
    x_range = [-6, 0];     % for doping or mu, depending on chosen mode

    % if true, ask folder by uigetdir; otherwise use default_root directly
    use_uigetdir = true;

    % ------------------------------------------------------------
    % choose x-axis mode
    % ------------------------------------------------------------
    choice = questdlg( ...
        'Choose x-axis variable:', ...
        'X-axis selection', ...
        'doping', 'mu', 'doping');

    if isempty(choice)
        error('No x-axis mode selected.');
    end

    switch lower(choice)
        case 'doping'
            x_mode = 'doping';
            x_label = 'Doping (10^{12} cm^{-2})';
        case 'mu'
            x_mode = 'mu';
            x_label = '\mu (eV)';
        otherwise
            error('Unknown x-axis selection.');
    end

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
    S = struct('name', {}, 'path', {}, 'D_meV', {}, ...
               'mu', {}, 'doping', {}, 'DOS', {});
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
        if isempty(M) || size(M,2) < 5
            fprintf('[skip] numeric block has fewer than 5 columns: %s\n', fname);
            continue;
        end

        % columns: i E(eV) filling doping DOS
        mu     = M(:,2);
        doping = M(:,4);
        DOS    = M(:,5);

        good = isfinite(mu) & isfinite(doping) & isfinite(DOS);
        mu     = mu(good);
        doping = doping(good);
        DOS    = DOS(good);

        if isempty(mu)
            fprintf('[skip] empty valid data: %s\n', fname);
            continue;
        end

        % choose x
        switch x_mode
            case 'doping'
                x = doping;
            case 'mu'
                x = mu;
        end

        % optional x-range restriction
        if use_x_range
            keep2 = (x >= x_range(1)) & (x <= x_range(2));
            mu     = mu(keep2);
            doping = doping(keep2);
            DOS    = DOS(keep2);
            x      = x(keep2);
        end

        if isempty(x)
            fprintf('[skip] empty after x-range cut: %s\n', fname);
            continue;
        end

        % sort by selected x for smooth curve
        [x_sorted, ord2] = sort(x);
        mu     = mu(ord2);
        doping = doping(ord2);
        DOS    = DOS(ord2);

        nkeep = nkeep + 1;
        S(nkeep).name   = fname;
        S(nkeep).path   = fpath;
        S(nkeep).D_meV  = D_meV;
        S(nkeep).mu     = mu;
        S(nkeep).doping = doping;
        S(nkeep).DOS    = DOS;
        S(nkeep).x      = x_sorted;
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

        if numel(pick_idx) < N_plot
            rest = setdiff(1:nall, pick_idx, 'stable');
            nneed = N_plot - numel(pick_idx);
            pick_idx = [pick_idx, rest(1:nneed)];
            pick_idx = sort(pick_idx);
        end
    end

    Splot = S(pick_idx);
    Dvals = [Splot.D_meV];

    fprintf('x-mode: %s\n', x_mode);
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

        plot(ax, Splot(i).x, Splot(i).DOS, ...
            'LineWidth', 1.5, ...
            'Color', thisColor, ...
            'DisplayName', sprintf('D = %.3f meV', Splot(i).D_meV));
    end

    xlabel(ax, x_label, 'FontSize', 16);
    ylabel(ax, 'DOS', 'FontSize', 16);
    title(ax, sprintf('DOS vs %s for D in [%g, %g] meV (%d curves)', ...
        x_mode, D_range(1), D_range(2), numel(Splot)), 'FontSize', 16);

    set(ax, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on');
    grid(ax, 'on');

    caxis(ax, [Dmin, Dmax]);
    cb = colorbar(ax);
    cb.Label.String = 'D (meV)';
    cb.FontSize = 14;

    % legend(ax, 'show', 'Location', 'best');

    % ------------------------------------------------------------
    % output
    % ------------------------------------------------------------
    out.root      = root;
    out.x_mode    = x_mode;
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