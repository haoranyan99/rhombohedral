function plot_hyst_metric_vs_artificial_polar()
clc; close all;

% =========================
% USER SETTINGS
% =========================
T_target = 6.5;          % manually set temperature
quantity = "psi";        % "psi" or "X"
FS = 16;

% =========================
% Select MAT file
% =========================
default_dir = 'E:\rg_master\rhombohedral\Matlab\hyst_data';

if ~isfolder(default_dir)
    default_dir = pwd;   % fallback
end

[fname, fpath] = uigetfile(fullfile(default_dir, '*.mat'), ...
    'Select hyst_allT_allArtificialPolar.mat');
if isequal(fname,0)
    return;
end

matfile = fullfile(fpath, fname);
S = load(matfile);

if ~isfield(S, "HYST")
    error("Selected MAT file does not contain variable HYST.");
end

HYST = S.HYST;

valid = arrayfun(@(x) isfield(x,"T") && isfinite(x.T) && ...
    isfield(x,"artificial_polar") && isfinite(x.artificial_polar), HYST);

entries = HYST(valid);

if isempty(entries)
    error("No valid HYST entries found.");
end

% =========================
% Select fixed temperature
% =========================
tolT = 1e-8;
E = entries(abs([entries.T] - T_target) < tolT);

if isempty(E)
    Tall = unique([entries.T]);
    error("No data found for T = %.6g K. Available T: %s", ...
        T_target, mat2str(Tall));
end

% =========================
% UI: choose plot type only
% =========================
plot_type_list = { ...
    'hyst strength: mean(abs(fwd - bwd))', ...
    'max value: max(abs fwd/bwd)'};

[idxM, ok] = listdlg( ...
    'PromptString', 'Select plot type:', ...
    'SelectionMode', 'single', ...
    'ListString', plot_type_list);

if ~ok
    return;
end

metric_id = idxM;

% =========================
% Compute metric
% =========================
polar = nan(numel(E),1);
metric = nan(numel(E),1);

for i = 1:numel(E)
    R = E(i);
    polar(i) = R.artificial_polar;

    [xf, yf, xb, yb] = get_quantity_curves_(R, quantity);

    okf = isfinite(xf) & isfinite(yf);
    okb = isfinite(xb) & isfinite(yb);

    xf = xf(okf);
    yf = yf(okf);
    xb = xb(okb);
    yb = yb(okb);

    if isempty(xf) || isempty(xb)
        continue;
    end

    switch metric_id
        case 1
            [xb, ord] = sort(xb);
            yb = yb(ord);

            [xb, ia] = unique(xb);
            yb = yb(ia);

            yb_on_fwd = interp1(xb, yb, xf, 'linear', NaN);

            ok2 = isfinite(yf) & isfinite(yb_on_fwd);

            if nnz(ok2) >= 3
                metric(i) = mean(abs(yf(ok2) - yb_on_fwd(ok2)));
            end

        case 2
            metric(i) = max([abs(yf(:)); abs(yb(:))], [], 'omitnan');
    end
end

ok = isfinite(polar) & isfinite(metric);
polar = polar(ok);
metric = metric(ok);

[polar, ord] = sort(polar);
metric = metric(ord);

if isempty(polar)
    error("No valid metric data found.");
end

% =========================
% Plot
% =========================
figure('Color','w','Position',[120 120 760 560]);

plot(polar, metric, 'o-', 'LineWidth', 2, 'MarkerSize', 7);
grid on; box on;

xlabel('Artificial polar / magnetic field proxy (meV)', 'Interpreter','latex');

if metric_id == 1
    ylabel(sprintf('Hysteresis strength of $|%s|$', quantity), 'Interpreter','latex');
    title_metric = 'hysteresis strength';
else
    ylabel(sprintf('Maximum value of $|%s|$', quantity), 'Interpreter','latex');
    title_metric = 'maximum value';
end

title(sprintf('$T = %.3f\\,\\mathrm{K}$, $|%s|$ %s vs artificial polar', ...
    T_target, quantity, title_metric), 'Interpreter','latex');

set(gca, ...
    'FontSize', FS, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth', 1.3, ...
    'TickDir','out');

fprintf("\nLoaded: %s\n", matfile);
fprintf("T = %.6g K\n", T_target);
fprintf("Quantity = %s\n", quantity);
fprintf("Metric = %s\n", title_metric);
disp(table(polar, metric));

end

function [xf, yf, xb, yb] = get_quantity_curves_(R, quantity)

xf = R.dop_fwd(:);
xb = R.dop_bwd(:);

switch lower(string(quantity))
    case "psi"
        if isfield(R,"psi_f_plot") && isfield(R,"psi_b_plot")
            yf = R.psi_f_plot(:);
            yb = R.psi_b_plot(:);
        elseif isfield(R,"abspsi_fwd") && isfield(R,"abspsi_bwd")
            yf = R.abspsi_fwd(:);
            yb = R.abspsi_bwd(:);
        elseif isfield(R,"psi_fwd") && isfield(R,"psi_bwd")
            yf = abs(R.psi_fwd(:));
            yb = abs(R.psi_bwd(:));
        else
            error("Cannot find psi data in HYST entry.");
        end

    case "x"
        if isfield(R,"X_f_plot") && isfield(R,"X_b_plot")
            yf = R.X_f_plot(:);
            yb = R.X_b_plot(:);
        elseif isfield(R,"absX_fwd") && isfield(R,"absX_bwd")
            yf = R.absX_fwd(:);
            yb = R.absX_bwd(:);
        elseif isfield(R,"X_fwd") && isfield(R,"X_bwd")
            yf = abs(R.X_fwd(:));
            yb = abs(R.X_bwd(:));
        else
            error("Cannot find X data in HYST entry.");
        end

    otherwise
        error("quantity must be psi or X.");
end

end