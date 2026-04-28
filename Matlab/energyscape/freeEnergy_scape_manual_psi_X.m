function plot_energy_landscape_psi_X_fixedLattice_withCurves()
% plot_energy_landscape_psi_X_fixedLattice_withCurves
%
% 2D free-energy landscape:
%   x-axis  : psi
%   y-axis  : X
%   colormap: hot
%
% Also plot two 1D cuts:
%   (1) F vs psi at fixed X = X_fix_for_psi_curve
%   (2) F vs X   at fixed psi = psi_fix_for_X_curve
%
% Electronic part:
%   F_el = 1/2 a1 psi^2 + b1/24 psi^4
%
% Fixed lattice part with three minima:
%   F_lat = kX * X^2 * (X^2 - X0^2)^2
%
% Coupling:
%   F_coup = lambda * psi * X

clc; close all;

%% =========================
% USER PARAMETERS
% ==========================

psi_fix_for_X_curve = -.564;
X_fix_for_psi_curve = 1.0272;

a1 = 0.1;
b1 = 12;
lambda = 0.4;

% psi_fix_for_X_curve = -1.015;
% X_fix_for_psi_curve = 1.0464;
% 
% a1 = -2;
% b1 = 12;
% lambda = 0.4;

% psi_fix_for_X_curve = 0;
% X_fix_for_psi_curve = 0;
% 
% a1 = 0.1;
% b1 = 12;
% lambda = 0;



% psi_fix_for_X_curve = -1.002;
% X_fix_for_psi_curve = 0;
% 
% a1 = -2;
% b1 = 12;
% lambda = 0;




% ---------------------------------
% manual cuts for 1D curves
% ---------------------------------



% fixed lattice shape
kX = 1;     % overall lattice energy scale
X0 = 1;     % side minima near X = +-X0

% plot range
psi_min = -1.5;
psi_max =  1.5;
X_min   = -1.2;
X_max   =  1.2;

% resolution
Npsi = 501;
NX   = 501;

% visualization
shift_min_to_zero   = true;
use_percentile_clip = true;
clip_percentile     = 98;
show_contour        = true;
num_contours        = 10;
FS = 16;

% local minima detection
exclude_boundary = true;   % ignore minima on outermost boundary
merge_tol_psi    = 0.02;   % merge nearby minima in psi
merge_tol_X      = 0.02;   % merge nearby minima in X


% whether to shift each 1D curve by its own minimum
shift_curve_min_to_zero = false;

%% =========================
% GRID
% ==========================
psiList = linspace(psi_min, psi_max, Npsi);
XList   = linspace(X_min, X_max, NX);
[PSI, XX] = meshgrid(psiList, XList);

dpsi = psiList(2) - psiList(1);
dX   = XList(2)   - XList(1);

%% =========================
% FREE ENERGY LANDSCAPE
% ==========================
F = free_energy_model(PSI, XX, a1, b1, lambda, kX, X0);

if shift_min_to_zero
    F = F - min(F(:));
end

Fplot = F;
if use_percentile_clip
    clip_val = prctile(Fplot(:), clip_percentile);
    Fplot(Fplot > clip_val) = clip_val;
end

%% =========================
% GLOBAL MINIMUM
% ==========================
[Fmin, idxMin] = min(F(:));
[rowMin, colMin] = ind2sub(size(F), idxMin);
psi_min_loc = PSI(rowMin, colMin);
X_min_loc   = XX(rowMin, colMin);

%% =========================
% FIND ALL LOCAL MINIMA ON GRID
% A point is a local minimum if it is <= all 8 neighbors
% and strictly < at least one neighbor.
% ==========================
isLocalMin = false(size(F));

row_start = 2;
row_end   = size(F,1)-1;
col_start = 2;
col_end   = size(F,2)-1;

if ~exclude_boundary
    row_start = 1;
    row_end   = size(F,1);
    col_start = 1;
    col_end   = size(F,2);
end

for i = row_start:row_end
    for j = col_start:col_end
        
        i1 = max(i-1, 1);
        i2 = min(i+1, size(F,1));
        j1 = max(j-1, 1);
        j2 = min(j+1, size(F,2));
        
        block = F(i1:i2, j1:j2);
        center = F(i,j);
        
        neigh = block(:);

        % center index inside the block
        nrow = i2 - i1 + 1;
        ncol = j2 - j1 + 1;
        center_idx = sub2ind([nrow, ncol], i - i1 + 1, j - j1 + 1);
        
        % remove center element
        neigh(center_idx) = [];
        
        if all(center <= neigh) && any(center < neigh)
            isLocalMin(i,j) = true;
        end
    end
end

psi_local = PSI(isLocalMin);
X_local   = XX(isLocalMin);
F_local   = F(isLocalMin);

%% =========================
% MERGE VERY CLOSE MINIMA
% ==========================
if ~isempty(psi_local)
    keep = true(size(psi_local));
    
    % sort by energy first so lower minima are kept preferentially
    [~, order] = sort(F_local, 'ascend');
    psi_local = psi_local(order);
    X_local   = X_local(order);
    F_local   = F_local(order);
    
    for n = 1:length(psi_local)
        if ~keep(n), continue; end
        
        for m = n+1:length(psi_local)
            if ~keep(m), continue; end
            
            if abs(psi_local(n) - psi_local(m)) < max(merge_tol_psi, 1.5*dpsi) && ...
               abs(X_local(n)   - X_local(m))   < max(merge_tol_X,   1.5*dX)
                keep(m) = false;
            end
        end
    end
    
    psi_local = psi_local(keep);
    X_local   = X_local(keep);
    F_local   = F_local(keep);
end

%% =========================
% PRINT LOCAL MINIMA
% ==========================
fprintf('\n==== Global minimum ====\n');
fprintf('psi = %+8.5f, X = %+8.5f, F = %+12.8f\n', ...
    psi_min_loc, X_min_loc, Fmin);

fprintf('\n==== Local minima found: %d ====\n', numel(F_local));
for n = 1:numel(F_local)
    fprintf('min %2d: psi = %+8.5f, X = %+8.5f, F = %+12.8f\n', ...
        n, psi_local(n), X_local(n), F_local(n));
end

%% =========================
% 1D CURVES
% ==========================
F_vs_psi = free_energy_model(psiList, X_fix_for_psi_curve, ...
    a1, b1, lambda, kX, X0);

F_vs_X = free_energy_model(psi_fix_for_X_curve, XList, ...
    a1, b1, lambda, kX, X0);

if shift_curve_min_to_zero
    F_vs_psi = F_vs_psi - min(F_vs_psi);
    F_vs_X   = F_vs_X   - min(F_vs_X);
end

% values of cut-lines on the 2D landscape
F_cut_psi = free_energy_model(psiList, X_fix_for_psi_curve, ...
    a1, b1, lambda, kX, X0);

F_cut_X = free_energy_model(psi_fix_for_X_curve, XList, ...
    a1, b1, lambda, kX, X0);

if shift_min_to_zero
    F_cut_psi = F_cut_psi - min(F(:));
    F_cut_X   = F_cut_X   - min(F(:));
end

%% =========================
% PLOT: 2D landscape
% ==========================
figure('Color','w','Position',[100 80 850 700]);

imagesc(psiList, XList, Fplot);
set(gca, 'YDir', 'normal');
axis tight;
colormap(hot);
colorbar;
hold on;

if show_contour
    contour(psiList, XList, F, num_contours, 'k-', 'LineWidth', 0.45);
end

% plot all local minima
plot(psi_local, X_local, 'bo', ...
    'MarkerSize', 7, 'LineWidth', 1.5, 'MarkerFaceColor', 'b');

% mark global minimum
plot(psi_min_loc, X_min_loc, 'co', ...
    'MarkerSize', 9, 'LineWidth', 1.8);

% draw the two cut lines
plot(psiList, X_fix_for_psi_curve * ones(size(psiList)), 'w--', 'LineWidth', 1.4);
plot(psi_fix_for_X_curve * ones(size(XList)), XList, 'c--', 'LineWidth', 1.4);

xlabel('\psi', 'FontSize', FS);
ylabel('X', 'FontSize', FS);
title('2D free-energy landscape', 'FontSize', FS);
set(gca, 'FontSize', FS, 'LineWidth', 1.2, 'Box', 'on');

%% =========================
% PLOT: F vs psi at fixed X
% ==========================
figure('Color','w','Position',[120 100 760 560]);
plot(psiList, F_vs_psi, 'k-', 'LineWidth', 2); hold on;

% mark local minima lying close to this fixed X cut
tol_X_cut = max(merge_tol_X, 1.5*dX);
mask_cut_psi = abs(X_local - X_fix_for_psi_curve) < tol_X_cut;
if any(mask_cut_psi)
    psi_cut_pts = psi_local(mask_cut_psi);
    F_cut_pts   = free_energy_model(psi_cut_pts, X_fix_for_psi_curve, ...
        a1, b1, lambda, kX, X0);
    if shift_curve_min_to_zero
        F_cut_pts = F_cut_pts - min(F_vs_psi);
    end
    plot(psi_cut_pts, F_cut_pts, 'bo', 'MarkerSize', 7, ...
        'LineWidth', 1.5, 'MarkerFaceColor', 'b');
end

xlabel('\psi', 'FontSize', FS);
ylabel('F(\psi)', 'FontSize', FS);
title(sprintf('1D cut: X = %.4f', X_fix_for_psi_curve), 'FontSize', FS);
set(gca, 'FontSize', FS, 'LineWidth', 1.2, 'Box', 'on');

%% =========================
% PLOT: F vs X at fixed psi
% ==========================
figure('Color','w','Position',[140 120 760 560]);
plot(XList, F_vs_X, 'k-', 'LineWidth', 2); hold on;

% mark local minima lying close to this fixed psi cut
tol_psi_cut = max(merge_tol_psi, 1.5*dpsi);
mask_cut_X = abs(psi_local - psi_fix_for_X_curve) < tol_psi_cut;
if any(mask_cut_X)
    X_cut_pts = X_local(mask_cut_X);
    F_cut_pts = free_energy_model(psi_fix_for_X_curve, X_cut_pts, ...
        a1, b1, lambda, kX, X0);
    if shift_curve_min_to_zero
        F_cut_pts = F_cut_pts - min(F_vs_X);
    end
    plot(X_cut_pts, F_cut_pts, 'bo', 'MarkerSize', 7, ...
        'LineWidth', 1.5, 'MarkerFaceColor', 'b');
end

xlabel('X', 'FontSize', FS);
ylabel('F(X)', 'FontSize', FS);
title(sprintf('1D cut: \\psi = %.4f', psi_fix_for_X_curve), 'FontSize', FS);
set(gca, 'FontSize', FS, 'LineWidth', 1.2, 'Box', 'on');

end

%% =========================
% LOCAL FUNCTION
% ==========================
function F = free_energy_model(psi, X, a1, b1, lambda, kX, X0)

F_el  = 0.5 .* a1 .* psi.^2 + (b1/24) .* psi.^4;
F_lat = kX .* X.^2 .* (X.^2 - X0.^2).^2;
F_cpl = lambda .* psi .* X;

F = F_el + F_lat + F_cpl;
end