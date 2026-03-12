function freeEnergy_scape_artificial()
    % 参数范围设置
    param_range = struct(...
        'T', [50, 120], ...
        'n', [-10, 10], ...
        'lambda', [0, 1]);

    % 自由能参数对 T,n 的依赖
    beta = @(T) 2 - (T-70)/20;
    n_1  = @(n) (n - abs(n))/2 + (abs(n) + n)/2; % == n

    a1_formula = @(T, n) -2.5 * beta(T) - 1.7*(n_1(n-2)+1.3) + 2;
    b1_formula = @(T) 1;
    b2_formula = @(T) -1 * beta(T);
    a2_formula = @(T,n) 0.5 * (beta(T).^2) + 1.4 - 0.075*(n_1(n-2)+1.3);
    c2_formula = @(T) 1;

    % 初始值
    psi1_init = 0;
    psi2_init = 0;

    % 创建 UI（加宽给右侧）
    fig = uifigure('Name', 'Interactive Free Energy', 'Position', [100 100 980 600]);

    % ===================== 左：2D 自由能图 =====================
    ax2d = uiaxes(fig, 'Position', [50 200 520 400], 'FontSize', 12);
    ax2d.TickLabelInterpreter = 'latex';
    xlabel(ax2d, '$\psi_{\rm CDW}$', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel(ax2d, '$\psi_{\rm lattice}$', 'FontSize', 12, 'Interpreter', 'latex');
    title(ax2d, 'Free Energy $F(\psi_1,\psi_2)$', 'FontSize', 16, 'Interpreter', 'latex');
    xlim(ax2d, [-12,12]);
    ylim(ax2d, [-7,7]);

    % ===================== 右：1D cuts（忽略耦合） =====================
    axPsi = uiaxes(fig, 'Position', [600 410 340 190], 'FontSize', 12);
    axPsi.TickLabelInterpreter = 'latex';
    xlabel(axPsi, '$\psi_1$', 'Interpreter', 'latex');
    ylabel(axPsi, '$F_\psi(\psi_1)-\min$', 'Interpreter', 'latex');
    title(axPsi, title_psi_(), 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'normal');

    axX = uiaxes(fig, 'Position', [600 200 340 190], 'FontSize', 12);
    axX.TickLabelInterpreter = 'latex';
    xlabel(axX, '$\psi_2$', 'Interpreter', 'latex');
    ylabel(axX, '$F_X(\psi_2)-\min$', 'Interpreter', 'latex');
    title(axX, title_X_(), 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'normal');

    % ---- 右侧系数显示框（放在右下，两个图下面）----
    coeff_box = uitextarea(fig, ...
        'Position', [600 40 340 140], ...
        'Editable', 'off', ...
        'FontSize', 12);

    % ===================== sliders（保留原风格） =====================
    y1 = 90;
    dy = 50;

    % T slider
    T_slider = uislider(fig, 'Position', [100 y1+2*dy 400 3], 'Limits', param_range.T, 'Value', 120);
    T_slider.MajorTicks = 0:10:120;
    T_slider.MinorTicks = [];

    % n slider
    n_slider = uislider(fig, 'Position', [100 y1+dy 400 3], 'Limits', param_range.n, 'Value', 1);
    n_slider.MajorTicks = -10:10;
    n_slider.MinorTicks = [];

    % lambda slider
    lambda_slider = uislider(fig, 'Position', [100 y1 400 3], 'Limits', param_range.lambda, 'Value', 0.5);
    lambda_slider.MajorTicks = 0:0.2:1;
    lambda_slider.MinorTicks = [];

    % labels
    uilabel(fig, 'Position', [70 y1+2*dy-10 50 22], 'Text', '$T$', 'FontSize', 12, 'Interpreter', 'latex');
    uilabel(fig, 'Position', [70 y1+dy-10 50 22], 'Text', '$n$', 'FontSize', 12, 'Interpreter', 'latex');
    uilabel(fig, 'Position', [70 y1-10 50 22], 'Text', '$\lambda$', 'FontSize', 12, 'Interpreter', 'latex');

    % ===================== 手动输入 + reset =====================
    psi1_input = uieditfield(fig, 'numeric', 'Position', [180 20 60 30], 'Value', psi1_init, 'FontSize', 12);
    psi2_input = uieditfield(fig, 'numeric', 'Position', [250 20 60 30], 'Value', psi2_init, 'FontSize', 12);
    uilabel(fig, 'Position', [180 50 60 22], 'Text', '$\psi_1$', 'FontSize', 12, 'Interpreter', 'latex');
    uilabel(fig, 'Position', [250 50 60 22], 'Text', '$\psi_2$', 'FontSize', 12, 'Interpreter', 'latex');

    reset_button = uibutton(fig, 'push', 'Position', [320 20 100 30], 'Text', 'Reset', 'FontSize', 12);

    % ===================== 全局最优点（延续迭代） =====================
    global psi1_opt psi2_opt;
    psi1_opt = psi1_init;
    psi2_opt = psi2_init;

    % slider 回调：参数变化 -> 重新极小化 -> 更新三幅图
    T_slider.ValueChangedFcn      = @(~,~) update_minimization();
    n_slider.ValueChangedFcn      = @(~,~) update_minimization();
    lambda_slider.ValueChangedFcn = @(~,~) update_minimization();

    % reset：用输入框作为初值继续最小化
    reset_button.ButtonPushedFcn = @(~,~) reset_initial_conditions();

    % 初始绘图
    update_plot();

    % ===================== nested functions =====================
    function reset_initial_conditions()
        psi1_opt = psi1_input.Value;
        psi2_opt = psi2_input.Value;
        update_minimization();
    end

    function update_minimization()
        T      = T_slider.Value;
        n      = n_slider.Value;
        lambda = lambda_slider.Value;

        a1 = a1_formula(T, n);
        b1 = b1_formula(T);
        a2 = a2_formula(T, n);
        b2 = b2_formula(T);
        c2 = c2_formula(T);

        [psi1_opt, psi2_opt, ~] = minimize_free_energy( ...
            a1, b1, a2, b2, c2, lambda, ...
            psi1_opt + 1e-3*rand(), psi2_opt + 1e-3*rand(), false);

        update_plot();
    end

    function update_plot()
        T      = T_slider.Value;
        n      = n_slider.Value;
        lambda = lambda_slider.Value;

        a1 = a1_formula(T, n);
        b1 = b1_formula(T);
        a2 = a2_formula(T, n);
        b2 = b2_formula(T);
        c2 = c2_formula(T);

        % -------- 2D (WITH coupling) --------
        F2D = free_energy(a1, b1, a2, b2, c2, lambda);

        psi1_vals = linspace(-12, 12, 120);
        psi2_vals = linspace(-7, 7, 120);
        [Psi1, Psi2] = meshgrid(psi1_vals, psi2_vals);
        F_vals = arrayfun(F2D, Psi1, Psi2);

        cla(ax2d);
        pcolor(ax2d, Psi1, Psi2, F_vals);
        shading(ax2d, 'interp');
        colormap(ax2d, 'turbo');
        colorbar(ax2d, 'Location', 'eastoutside', 'FontSize', 12);

        hold(ax2d, 'on');
        contour(ax2d, Psi1, Psi2, F_vals, 20, 'LineColor', 'k', 'LineWidth', 0.5);
        plot(ax2d, psi1_opt, psi2_opt, 'ro', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
        hold(ax2d, 'off');

        title(ax2d, sprintf('Free Energy (T=%.2f, n=%.2f, \\lambda=%.2f)', T, n, lambda), 'FontSize', 14);

        % -------- 1D cuts (IGNORE coupling) --------
        Fpsi = @(psi) 0.5*a1*psi.^2 + (1/factorial(4))*b1*psi.^4;
        FX   = @(X)   0.5*a2*X.^2   + (1/factorial(4))*b2*X.^4 + (1/factorial(6))*c2*X.^6;

        psi_line = linspace(-12, 12, 801);
        X_line   = linspace(-7, 7, 801);

        Fpsi_v = Fpsi(psi_line);
        FX_v   = FX(X_line);

        % shift by min (保持你现在 y label 的含义)
        Fpsi_min = min(Fpsi_v);
        FX_min   = min(FX_v);

        Fpsi_plot = Fpsi_v - Fpsi_min;
        FX_plot   = FX_v   - FX_min;

        % plot Fpsi
        cla(axPsi);
        plot(axPsi, psi_line, Fpsi_plot, 'LineWidth', 2);
        grid(axPsi, 'on'); box(axPsi, 'on');
        hold(axPsi, 'on');
        plot(axPsi, psi1_opt, Fpsi(psi1_opt) - Fpsi_min, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        hold(axPsi, 'off');
        title(axPsi, title_psi_(), 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'normal');

        % plot FX
        cla(axX);
        plot(axX, X_line, FX_plot, 'LineWidth', 2);
        grid(axX, 'on'); box(axX, 'on');
        hold(axX, 'on');
        plot(axX, psi2_opt, FX(psi2_opt) - FX_min, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        hold(axX, 'off');
        title(axX, title_X_(), 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'normal');

        % -------- update coefficient box --------
        coeff_box.Value = { ...
            sprintf('T = %.6g,   n = %.6g,   lambda = %.6g', T, n, lambda), ...
            ' ', ...
            sprintf('a1 = %.6g', a1), ...
            sprintf('b1 = %.6g', b1), ...
            sprintf('a2 = %.6g', a2), ...
            sprintf('b2 = %.6g', b2), ...
            sprintf('c2 = %.6g', c2), ...
            ' ', ...
            sprintf('psi1* = %.6g,   psi2* = %.6g', psi1_opt, psi2_opt) ...
        };
    end
end

% ---------- titles (as you requested) ----------
function t = title_psi_()
    t = sprintf('$F_\\psi=\\frac12 a_1\\psi^2+\\frac{1}{4!}b_1\\psi^4\\quad$');
end

function t = title_X_()
    t = sprintf('$F_X=\\frac12 a_2X^2+\\frac{1}{4!}b_2X^4+\\frac{1}{6!}c_2X^6\\quad$');
end
