function C = detect_chi_columns(M, fpath, iq_pick, jq_pick)
%DETECT_CHI_COLUMNS Return iq/jq/Re(chi) columns for chi txt tables.
% Supports the current format:
%   # iq jq qx qy chi_real chi_imag nKpair nK
% and older tables with a leading idx column.

    C = detect_from_header_(fpath);
    if ~isempty(C)
        return;
    end

    candidates = [
        struct('iq',1,'jq',2,'Re',5)
        struct('iq',2,'jq',3,'Re',6)
    ];

    for i = 1:numel(candidates)
        cand = candidates(i);
        if size(M,2) < max([cand.iq, cand.jq, cand.Re])
            continue;
        end

        if nargin >= 4 && any(M(:,cand.iq)==iq_pick & M(:,cand.jq)==jq_pick)
            C = cand;
            return;
        end
    end

    if size(M,2) >= 6
        C = struct('iq',1,'jq',2,'Re',5);
    else
        error("detect_chi_columns: unsupported chi table with %d columns: %s", ...
            size(M,2), fpath);
    end
end

function C = detect_from_header_(fpath)

    C = [];

    fid = fopen(fpath,'r');
    if fid < 0
        return;
    end

    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    for t = 1:3000
        ln = fgetl(fid);
        if ~ischar(ln)
            break;
        end

        s = strtrim(string(ln));
        if ~startsWith(s, "#")
            break;
        end

        s = strtrim(erase(s, "#"));
        slo = lower(s);

        if ~contains(slo, "iq") || ~contains(slo, "jq") || ...
           ~(contains(slo, "chi_real") || contains(slo, "chi_re"))
            continue;
        end

        toks = regexp(char(s), "\s+", "split");
        toks = toks(~cellfun(@isempty, toks));
        if isempty(toks)
            continue;
        end

        toks_norm = lower(string(toks));
        toks_norm = regexprep(toks_norm, "[^a-z0-9_]", "");
        if ~isempty(toks_norm) && toks_norm(1) == "columns"
            toks_norm = toks_norm(2:end);
        end

        iq_col = find(toks_norm == "iq", 1);
        jq_col = find(toks_norm == "jq", 1);

        re_names = ["chi_real", "chi_re", "chire", "real", "re"];
        re_col = [];
        for name = re_names
            re_col = find(toks_norm == name, 1);
            if ~isempty(re_col)
                break;
            end
        end

        if ~isempty(iq_col) && ~isempty(jq_col) && ~isempty(re_col)
            C = struct('iq',iq_col,'jq',jq_col,'Re',re_col);
            return;
        end
    end
end
