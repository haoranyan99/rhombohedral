function S = plot_hyst_width_from_one_mat()
clc; close all;

sep_threshold = 1e-2;

[fname, fpath] = uigetfile("*.mat", "Select hyst_allT.mat");
if isequal(fname,0), error("User cancelled"); end

Sdata = load(fullfile(fpath, fname));
HYST = Sdata.HYST;

nT = numel(HYST);

T_list = nan(nT,1);
Width = nan(nT,1);
Left  = nan(nT,1);
Right = nan(nT,1);

for i = 1:nT

    H = HYST(i);

    T_list(i) = H.T;

    dop_fwd = H.dop_fwd(:);
    dop_bwd = H.dop_bwd(:);

    psi_f = H.psi_f_plot(:);
    psi_b = H.psi_b_plot(:);

    [dop_b_sort, ord] = sort(dop_bwd);
    psi_b_sort = psi_b(ord);

    psi_b_interp = interp1(dop_b_sort, psi_b_sort, dop_fwd, ...
        'linear','extrap');

    sep = abs(psi_f - psi_b_interp);

    idx = find(sep > sep_threshold);

    if isempty(idx)
        Width(i) = 0;
    else
        Left(i)  = dop_fwd(idx(1));
        Right(i) = dop_fwd(idx(end));
        Width(i) = abs(Right(i) - Left(i));
    end

end

% ===== sort =====
[T_list, ord] = sort(T_list);
Width = Width(ord);

% ===== plot =====
figure('Color','w','Position',[120 120 760 560]);
plot(T_list, Width, 'o-', 'LineWidth',2);

grid on; box on;

xlabel('$T$ (K)','Interpreter','latex');
ylabel('$\Delta n_{\rm hyst}$','Interpreter','latex');

title(sprintf('Hysteresis width (threshold=%.1e)', sep_threshold), ...
    'Interpreter','latex');

set(gca,'FontSize',16,'TickLabelInterpreter','latex');

S = struct();
S.T = T_list;
S.Width = Width;

end