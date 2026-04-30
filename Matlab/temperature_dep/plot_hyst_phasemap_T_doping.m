function plot_hyst_phasemap_T_doping()
clc; close all;

% =========================
% USER SETTINGS
% =========================
polar_target = 0;       % 固定磁场 (meV)
quantity = "psi";       % "psi" or "X"
FS = 16;

doping_N = 300;
use_abs_diff = true;

default_dir = 'E:\rg_master\rhombohedral\Matlab\hyst_data';

% =========================
% Load MAT
% =========================
if ~isfolder(default_dir)
    default_dir = pwd;
end

[fname, fpath] = uigetfile(fullfile(default_dir, '*.mat'), ...
    'Select hyst_allT_allArtificialPolar.mat');

if isequal(fname,0); return; end

S = load(fullfile(fpath, fname));
HYST = S.HYST;

valid = arrayfun(@(x) isfield(x,"T") && isfield(x,"artificial_polar"), HYST);
entries = HYST(valid);

% =========================
% 固定 polar
% =========================
tolP = 1e-8;
E = entries(abs([entries.artificial_polar] - polar_target) < tolP);

if isempty(E)
    Pall = unique([entries.artificial_polar]);
    error("No data for polar=%.6g. Available: %s", polar_target, mat2str(Pall));
end

% =========================
% 统一 doping grid
% =========================
dop_min = -inf;
dop_max = inf;

for i = 1:numel(E)
    [xf,~,xb,~] = get_quantity(E(i), quantity);

    xf = xf(isfinite(xf));
    xb = xb(isfinite(xb));

    if isempty(xf) || isempty(xb); continue; end

    dop_min = max(dop_min, max(min(xf), min(xb)));
    dop_max = min(dop_max, min(max(xf), max(xb)));
end

dop_grid = linspace(dop_min, dop_max, doping_N);

% =========================
% 计算 phase map
% =========================
Tlist = nan(numel(E),1);
Z = nan(numel(E), doping_N);

for i = 1:numel(E)
    R = E(i);
    Tlist(i) = R.T;

    [xf,yf,xb,yb] = get_quantity(R, quantity);

    okf = isfinite(xf)&isfinite(yf);
    okb = isfinite(xb)&isfinite(yb);

    xf = xf(okf); yf = yf(okf);
    xb = xb(okb); yb = yb(okb);

    if numel(xf)<3 || numel(xb)<3; continue; end

    [xf,ia]=unique(xf); yf=yf(ia);
    [xb,ia]=unique(xb); yb=yb(ia);

    yf_grid = interp1(xf,yf,dop_grid,'linear',NaN);
    yb_grid = interp1(xb,yb,dop_grid,'linear',NaN);

    if use_abs_diff
        Z(i,:) = abs(yf_grid - yb_grid);
    else
        Z(i,:) = yf_grid - yb_grid;
    end
end

% 排序
[Tlist,ord] = sort(Tlist);
Z = Z(ord,:);

% =========================
% Plot
% =========================
figure('Color','w','Position',[120 120 850 620]);

imagesc(dop_grid, Tlist, Z);
set(gca,'YDir','normal');

colormap(flipud(hot));
caxis([0, max(Z(:),[],'omitnan')]);

cb = colorbar;
cb.Label.String = sprintf('$||%s|_{fwd}-|%s|_{bwd}|$',quantity,quantity);
cb.Label.Interpreter = 'latex';

xlabel('Doping $(10^{12}\,\mathrm{cm}^{-2})$','Interpreter','latex');
ylabel('Temperature (K)','Interpreter','latex');

title(sprintf('Hysteresis phase map (%s), polar=%.3f meV', ...
    quantity, polar_target), 'Interpreter','latex');

set(gca,'FontSize',FS,'TickLabelInterpreter','latex','LineWidth',1.3);

end

% =========================
% helper
% =========================
function [xf,yf,xb,yb] = get_quantity(R, quantity)

xf = R.dop_fwd(:);
xb = R.dop_bwd(:);

switch lower(string(quantity))
    case "psi"
        yf = R.psi_f_plot(:);
        yb = R.psi_b_plot(:);
    case "x"
        yf = R.X_f_plot(:);
        yb = R.X_b_plot(:);
end

end