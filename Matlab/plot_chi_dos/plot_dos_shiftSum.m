function out_png_list = plot_dos_shiftSum()
% UI choose one or multiple DOS files
% Plot shift-sum DOS vs E
% NEW: click any point -> datatip shows doping (from DATA column)

out_png_list = strings(0);

% ==========================
% USER SETTINGS
% ==========================
start_dir = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\";
mu_shift_meV = 1;     % meV 
overlay_original = true;    
save_png = false;        % optional
FS = 14;

% ==========================
% Select files
% ==========================
[file, path] = uigetfile("*.txt", "Select DOS txt file(s)", start_dir, "MultiSelect","on");
if isequal(file,0), error("No file selected."); end
if ischar(file), file = {file}; end

% ==========================
% Loop
% ==========================
for k = 1:numel(file)
    in_path = fullfile(path, file{k});
    [~, base, ~] = fileparts(in_path);

    raw = readmatrix(in_path, "FileType","text", "CommentStyle","#");
    if isempty(raw) || size(raw,2) < 3
        fprintf("[skip] bad file: %s\n", in_path);
        continue;
    end

    % ---- get E, DOS, doping from DATA columns ----
    % Typical: i  E  filling  doping  DOS  (>=5 cols)
    ncol = size(raw,2);
    E = raw(:,2);

    if ncol >= 5
        doping = raw(:,4);
        DOS = raw(:,5);
    elseif ncol == 4
        doping = raw(:,3);
        DOS = raw(:,4);
    else
        fprintf("[skip] no doping column: %s\n", in_path);
        continue;
    end

    ok = isfinite(E) & isfinite(DOS) & isfinite(doping);
    E = E(ok); DOS = DOS(ok); doping = doping(ok);
    if numel(E) < 3
        fprintf("[skip] too few points: %s\n", in_path);
        continue;
    end

    % sort by E (keep doping aligned!)
    [Egrid, ord] = sort(E);
    DOSgrid = DOS(ord);
    dopgrid = doping(ord);

    % shift-sum (keep SAME Egrid; doping for datatip still uses dopgrid)
    shift_eV = 0.001 * mu_shift_meV;
    DOS_left  = interp1(Egrid, DOSgrid, Egrid + shift_eV, 'linear', NaN);
    DOS_right = interp1(Egrid, DOSgrid, Egrid - shift_eV, 'linear', NaN);
    DOS_sum   = 0.5*(DOS_left + DOS_right);

    % ==========================
    % Plot
    % ==========================
    fig = figure('Color','w','Units','pixels','Position',[120 120 820 560], ...
        'Name', "shift-sum DOS (click to show doping)");

    ax = axes(fig); hold(ax,'on');

    % shift-sum curve
    hSum = plot(ax, Egrid, DOS_sum, '--', 'LineWidth', 2.2, ...
        'DisplayName', sprintf("shift-sum, \\mu=%.4g meV", mu_shift_meV));

    % store doping for datatip
    setappdata(hSum, "dopgrid", dopgrid);
    setappdata(hSum, "Egrid",   Egrid);

    if overlay_original
        hOrg = plot(ax, Egrid, DOSgrid, 'LineWidth', 2.0, 'DisplayName', "original DOS");
        setappdata(hOrg, "dopgrid", dopgrid);
        setappdata(hOrg, "Egrid",   Egrid);
        legend(ax,'Location','best');
    end

    xlabel(ax,'E (eV)','FontSize',FS);
    ylabel(ax,'DOS','FontSize',FS);
    title(ax, base + "  (click point to show doping)", 'Interpreter','none');
    set(ax,'FontSize',FS,'LineWidth',1,'Box','on','TickDir','out');
    grid(ax,'on');

    % enable datacursor + custom tip
    dcm = datacursormode(fig);
    set(dcm, "Enable","on", "DisplayStyle","window", "SnapToDataVertex","on");
    set(dcm, "UpdateFcn", @tip_show_doping);

    if save_png
        out_png = fullfile(path, base + "_shiftSum_clickDoping.png");
        exportgraphics(fig, out_png,'Resolution',300);
        out_png_list(end+1) = string(out_png);
    end
end

end

% ==========================
% DataTip callback
% ==========================
function txt = tip_show_doping(~, event)
    % event.Target is the line handle
    h = event.Target;

    % Base info from clicked point
    pos = event.Position; % [E, DOS, (maybe z)]
    Eclick = pos(1);
    DOSclick = pos(2);

    % Default
    txt = {
        sprintf("E = %.10g eV", Eclick)
        sprintf("DOS = %.10g",  DOSclick)
    };

    % Try to fetch aligned doping via DataIndex (best)
    dopgrid = [];
    try
        dopgrid = getappdata(h, "dopgrid");
    catch
    end

    if ~isempty(dopgrid)
        idx = [];
        try
            idx = event.DataIndex; % MATLAB provides this for line/scatter
        catch
        end

        if ~isempty(idx) && idx >= 1 && idx <= numel(dopgrid)
            txt{end+1} = sprintf("doping = %.10g", dopgrid(idx));
            return;
        end
    end

    % Fallback: nearest by E
    Egrid = [];
    try  
        Egrid = getappdata(h, "Egrid");
    catch
    end  
    if ~isempty(Egrid) && ~isempty(dopgrid) && numel(Egrid)==numel(dopgrid)
        [~, ii] = min(abs(Egrid - Eclick));
        txt{end+1} = sprintf("doping = %.10g", dopgrid(ii));
    end
end