function out_png = plot_bareband()

    % -----------------------------
    % choose file
    % -----------------------------
    default_root = "/Users/haoranyan/rg_master/data/";
    if ~isfolder(default_root)
        warning("Default root does not exist: %s\nFallback to current folder.", default_root);
        default_root = pwd;
    end

    default_path = fullfile(default_root, "*.txt");
    [fname, fpath] = uigetfile({'*.txt','DOS data (*.txt)'}, ...
                               'Select RG DOS file', default_path);
    if isequal(fname,0)
        error("No file selected.");
    end

    in_path = fullfile(fpath, fname);
    if ~isfile(in_path)
        error("Input file not found: %s", in_path);
    end

    % -----------------------------
    % read numeric data
    % -----------------------------
    A = readmatrix(in_path, 'CommentStyle', '#');
    if isempty(A) || size(A,2) < 5
        error('Data seems empty or has too few columns (need at least 5).');
    end

    s  = A(:,2);
    E  = A(:, 5:end);
    norb = size(E,2);
    fprintf("Detected %d bands.\n", norb);

    % -----------------------------
    % parse xticks from header
    % -----------------------------
    xtick_pos = [];
    xtick_lab = {};

    fid = fopen(in_path, 'r');
    if fid < 0
        error('Cannot open %s', in_path);
    end

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end
        line = strtrim(line);

        if startsWith(line, '# xticks:')
            blocks = regexp(line, '\(([^)]+)\)', 'tokens');
            for i = 1:numel(blocks)
                item = strtrim(blocks{i}{1});   % e.g. "K,0.123"
                parts = strsplit(item, ',');
                if numel(parts) >= 2
                    lab = strtrim(parts{1});
                    pos = str2double(strtrim(parts{2}));
                    if ~isnan(pos)
                        xtick_lab{end+1} = lab; %#ok<AGROW>
                        xtick_pos(end+1) = pos; %#ok<AGROW>
                    end
                end
            end
            break;
        end
    end
    fclose(fid);

    % fallback if header missing
    if isempty(xtick_pos)
        xtick_pos = linspace(s(1), s(end), 3);
        xtick_lab = {'start','mid','end'};
    end

    % -----------------------------
    % plot
    % -----------------------------
    fig = figure('Color','w');
    set(fig, 'Units','pixels', 'Position',[120,120,900,520]);

    hold on; box on; grid on;
    for ib = 1:norb
        plot(s, E(:,ib), 'LineWidth', 1.1);
    end

    xlabel('k-line coordinate s (1/\AA)');
    ylabel('Energy (eV)');
    set(gca, 'FontSize', 12);
    set(gca, 'XTick', xtick_pos, 'XTickLabel', xtick_lab);

    % -----------------------------
    % save
    % -----------------------------
    if ~exist('plot','dir'); mkdir('plot'); end
    [~, base, ~] = fileparts(in_path);
    out_png = fullfile('plot', base + ".png");
    saveas(fig, out_png);

    fprintf("Saved figure: %s\n", out_png);
end
