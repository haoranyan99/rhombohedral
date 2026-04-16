function plot_energy_landscape_psi_X_fixedLattice()
% plot_energy_landscape_psi_X_fixedLattice
% 2D free-energy landscape:
% x-axis  : psi
% y-axis  : X
% colormap: hot
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
% a1 = -2;
% b1 = 12;
% lambda = 0.2;

a1 = 0.1;
b1 = 12;
lambda = 0.2;

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

%% =========================
% GRID
% ==========================
psiList = linspace(psi_min, psi_max, Npsi);
XList   = linspace(X_min, X_max, NX);
[PSI, XX] = meshgrid(psiList, XList);

dpsi = psiList(2) - psiList(1);
dX   = XList(2)   - XList(1);

%% =========================
% FREE ENERGY
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
        center_idx = ceil(numel(neigh)/2);
        
        % remove center element
        neigh(center_idx) = [];
        
        if all(center <= neigh) && any(center < neigh)
            isLocalMin(i,j) = true;
        end
    end
end

[minRows, minCols] = find(isLocalMin);

psi_local = PSI(isLocalMin);
X_local   = XX(isLocalMin);
F_local   = F(isLocalMin);

%% =========================
% MERGE VERY CLOSE MINIMA
% (avoid many nearly identical minima due to dense mesh)
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
fprintf('\n==== Local minima found: %d ====\n', numel(F_local));
for n = 1:numel(F_local)
    fprintf('min %2d: psi = %+8.5f, X = %+8.5f, F = %+12.8f\n', ...
        n, psi_local(n), X_local(n), F_local(n));
end

%% =========================
% PLOT
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

% plot all local minima (same style)
plot(psi_local, X_local, 'bo', ...
    'MarkerSize', 7, 'LineWidth', 1.5, 'MarkerFaceColor', 'b');

xlabel('\psi', 'FontSize', FS);
ylabel('X', 'FontSize', FS);
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