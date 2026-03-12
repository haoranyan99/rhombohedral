function db = build_fermi_mat_database()
% BUILD_FERMI_MAT_DATABASE
% UI version for fermi data:
%   1) choose fermi root folder by popup
%   2) parse root folder name, e.g. fermi_sk_mu_400
%   3) recursively scan:
%        D{XXX}/T{XXX}/mu{XXXX}/fermiPatch*.txt
%   4) AFTER processing, choose output SAVE folder by popup
%   5) save everything into one .mat database
%
% Output:
%   db : returned database struct
%
% Example folder name:
%   fermi_sk_mu_400
% -> prefix = "fermi", model = "sk", mode = "mu", Nk = 400

    % ---------------- default folder ----------------
    default_folder = "/Users/haoranyan/data/rhombohedral";
    if ~isfolder(default_folder)
        default_folder = pwd;
    end

    % ---------------- choose root folder ----------------
    root_folder = uigetdir(default_folder, ...
        'Select fermi root folder to parse (e.g. fermi_sk_mu_400)');
    if isequal(root_folder, 0)
        fprintf('User canceled root-folder selection.\n');
        db = [];
        return;
    end
    root_folder = string(root_folder);

    % ---------------- parse root folder name ----------------
    [~, root_name] = fileparts(char(root_folder));
    root_name = string(root_name);

    meta_root = parse_fermi_root_folder_name(root_name);

    fprintf('Selected root folder:\n  %s\n', root_folder);
    fprintf('Parsed root metadata: prefix=%s, model=%s, mode=%s, Nk=%g\n', ...
        meta_root.prefix, meta_root.model, meta_root.mode, meta_root.Nk);

    % ---------------- recursively find fermiPatch*.txt ----------------
    files = dir(fullfile(root_folder, "**", "fermiPatch*.txt"));
    fprintf("Found %d fermiPatch*.txt files under:\n  %s\n", numel(files), root_folder);

    if isempty(files)
        warning('No fermiPatch*.txt files found.');
        db = [];
        return;
    end

    % ---------------- record template ----------------
    template = struct( ...
        'relpath', "", ...
        'folder', "", ...
        'name', "", ...
        'header', "", ...
        'data', [], ...
        'bytes', NaN, ...
        'modified', "", ...
        'D', NaN, ...
        'T', NaN, ...
        'mu', NaN, ...
        'model', "", ...
        'mode', "", ...
        'Nk', NaN ...
    );

    recs = repmat(template, numel(files), 1);

    % ---------------- waitbar ----------------
    wb = waitbar(0, 'Building FERMI MAT database...', 'Name', 'build_fermi_mat_database');
    cleanupObj = onCleanup(@() safe_close_waitbar(wb)); %#ok<NASGU>

    % ---------------- main loop ----------------
    for k = 1:numel(files)
        fullpath = string(fullfile(files(k).folder, files(k).name));
        relpath  = string(erase(fullpath, root_folder + filesep));

        [headerText, numericData] = read_txt_with_header(fullpath);
        [Dval, Tval, muVal] = parse_fermi_tree_path(relpath);

        recs(k).relpath   = relpath;
        recs(k).folder    = string(fileparts(char(relpath)));
        recs(k).name      = string(files(k).name);
        recs(k).header    = headerText;
        recs(k).data      = numericData;
        recs(k).bytes     = files(k).bytes;
        recs(k).modified  = string(files(k).date);

        recs(k).D         = Dval;
        recs(k).T         = Tval;
        recs(k).mu        = muVal;

        recs(k).model     = meta_root.model;
        recs(k).mode      = meta_root.mode;
        recs(k).Nk        = meta_root.Nk;

        if isgraphics(wb)
            waitbar(k / numel(files), wb, sprintf('Processing %d / %d', k, numel(files)));
        end

        if mod(k,100) == 0 || k == numel(files)
            fprintf("Processed %d / %d\n", k, numel(files));
        end
    end

    % ---------------- build db ----------------
    db = struct();
    db.root         = root_folder;
    db.root_name    = root_name;
    db.created_at   = string(datetime("now"));
    db.pattern      = "{root}/D{XXX}/T{XXX}/mu{XXXX}/fermiPatch*.txt";

    db.prefix       = meta_root.prefix;
    db.model        = meta_root.model;
    db.mode         = meta_root.mode;
    db.Nk           = meta_root.Nk;

    db.nfiles       = numel(recs);
    db.files        = recs;

    % ---------------- compact index ----------------
    db.index = struct();
    db.index.relpath = string({recs.relpath})';
    db.index.name    = string({recs.name})';
    db.index.D       = [recs.D]';
    db.index.T       = [recs.T]';
    db.index.mu      = [recs.mu]';
    db.index.model   = string({recs.model})';
    db.index.mode    = string({recs.mode})';
    db.index.Nk      = [recs.Nk]';

    % ---------------- optional table ----------------
    db.table = table( ...
        string({recs.relpath})', ...
        string({recs.name})', ...
        [recs.D]', ...
        [recs.T]', ...
        [recs.mu]', ...
        string({recs.model})', ...
        string({recs.mode})', ...
        [recs.Nk]', ...
        'VariableNames', {'relpath','name','D','T','mu','model','mode','Nk'});

    % ---------------- choose SAVE folder AFTER processing ----------------
    save_folder = uigetdir(default_folder, ...
        'Select folder to save FERMI MAT database');
    if isequal(save_folder, 0)
        fprintf('User canceled save-folder selection. Database kept in workspace only.\n');
        return;
    end
    save_folder = string(save_folder);

    out_mat = string(fullfile(save_folder, root_name + ".mat"));
    db.save_path = out_mat;

    % ---------------- overwrite protection ----------------
    if isfile(out_mat)
        choice = questdlg( ...
            sprintf('File already exists:\n%s\n\nOverwrite?', out_mat), ...
            'Overwrite existing MAT file?', ...
            'Overwrite', 'Cancel', 'Cancel');

        if ~strcmp(choice, 'Overwrite')
            fprintf('User canceled overwrite. Database kept in workspace only.\n');
            return;
        end
    end

    % ---------------- save ----------------
    save(out_mat, "db", "-v7.3");
    fprintf("Saved FERMI MAT database:\n  %s\n", out_mat);

    msgbox(sprintf('FERMI MAT database saved successfully:\n%s', out_mat), ...
           'Success', 'help');
end

function meta = parse_fermi_root_folder_name(root_name)
% Parse root folder name like:
%   fermi_sk_mu_400
%
% Output:
%   meta.prefix = "fermi"
%   meta.model  = "sk"
%   meta.mode   = "mu"
%   meta.Nk     = 400

    meta = struct();
    meta.prefix = "";
    meta.model  = "";
    meta.mode   = "";
    meta.Nk     = NaN;

    parts = split(string(root_name), "_");

    if numel(parts) >= 4
        meta.prefix = parts(1);
        meta.model  = parts(2);
        meta.mode   = parts(3);
        meta.Nk     = str2double(parts(4));
    else
        tok = regexp(char(root_name), ...
            '^([A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)_([-+]?\d+)$', ...
            'tokens', 'once');

        if ~isempty(tok)
            meta.prefix = string(tok{1});
            meta.model  = string(tok{2});
            meta.mode   = string(tok{3});
            meta.Nk     = str2double(tok{4});
        else
            warning('Root folder name "%s" does not match expected pattern fermi_sk_mu_400.', root_name);
        end
    end
end

function [headerText, numericData] = read_txt_with_header(filename)
% Read header lines (# or %) and numeric part.

    fid = fopen(filename, 'r');
    if fid < 0
        error("Cannot open file: %s", filename);
    end

    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    headerLines = {};
    numericLines = {};

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line)
            continue;
        end

        s = strtrim(line);
        if isempty(s) || startsWith(s, "#") || startsWith(s, "%")
            headerLines{end+1} = line; %#ok<AGROW>
        else
            numericLines{end+1} = line; %#ok<AGROW>
        end
    end

    headerText = string(strjoin(headerLines, newline));

    if isempty(numericLines)
        numericData = [];
        return;
    end

    tmp = [tempname, '.txt'];
    fid2 = fopen(tmp, 'w');
    if fid2 < 0
        error("Cannot create temporary file for numeric parsing.");
    end

    for i = 1:numel(numericLines)
        fprintf(fid2, '%s\n', numericLines{i});
    end
    fclose(fid2);

    numericData = readmatrix(tmp);
    delete(tmp);
end

function [Dval, Tval, muVal] = parse_fermi_tree_path(relpath)
% Parse path like:
%   D0.067/T4.000/mu1.005600/fermiPatch_xxx.txt

    relpath = char(relpath);

    Dval  = parse_one(relpath, '[/\\]D([-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?)');
    Tval  = parse_one(relpath, '[/\\]T([-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?)');
    muVal = parse_one(relpath, '[/\\]mu([-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?)');
end

function val = parse_one(str, pattern)
    tok = regexp(str, pattern, 'tokens', 'once');
    if isempty(tok)
        val = NaN;
    else
        val = str2double(tok{1});
    end
end

function safe_close_waitbar(wb)
    if ~isempty(wb) && isgraphics(wb)
        close(wb);
    end
end