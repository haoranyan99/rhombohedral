function [psi1_opt, psi2_opt, F_opt] = minimize_free_energy(a1, b1, a2, b2, c2, lambda, psi1_init, psi2_init, is_plot)
    tic;
    % 参数设置
    tol = 1e-6; % 迭代的容差
    max_iter = 1000; % 最大迭代次数
    iter_points = []; % 用于存储迭代过程中 (psi1, psi2) 点的轨迹
    total_secant_iter = 0; % 记录总的弦截法迭代次数
    cg_iter_count = 0; % 记录共轭梯度的迭代次数
    
    % 自由能
    [F, grad_F] = free_energy(a1, b1, a2, b2, c2, lambda);

    % 初始化共轭梯度法
    psi1 = psi1_init;
    psi2 = psi2_init;
    grad = grad_F(psi1, psi2);
    d = -grad;  % 初始搜索方向是负梯度
    iter = 0;

    % 存储初始点
    iter_points = [iter_points; psi1, psi2];

    % 迭代寻找最小值
    while norm(grad) > tol && iter < max_iter
        % 线搜索找到最优步长 alpha 并记录弦截法的迭代次数
        [alpha, secant_iter] = line_search(F, grad_F, psi1, psi2, d);
        total_secant_iter = total_secant_iter + secant_iter; % 累加弦截法的迭代次数

        % 更新 psi1 和 psi2
        psi1_new = psi1 + alpha * d(1);
        psi2_new = psi2 + alpha * d(2);

        % 计算新的梯度
        grad_new = grad_F(psi1_new, psi2_new);
        
        % Polak-Ribière+ CG
        beta = max((grad_new' * (grad_new - grad)) / (grad' * grad),0); 
        d_new = -grad_new + beta * d;
        if(dot(d_new,-grad_new)<0)
            d_new = -grad_new;
        end

        % 更新 psi 和梯度
        psi1 = psi1_new;
        psi2 = psi2_new;
        grad = grad_new;
        d = d_new;

        % 存储当前点
        iter_points = [iter_points; psi1, psi2];

        % 迭代计数
        iter = iter + 1;
    end
    
    % 记录共轭梯度的总迭代次数
    cg_iter_count = iter;

    % 输出结果
    psi1_opt = psi1;
    psi2_opt = psi2;
    F_opt = F(psi1_opt, psi2_opt);

    % 计算平均的弦截法迭代步数
    avg_secant_iter = total_secant_iter / cg_iter_count;

    % 绘制图像
    if is_plot
        fprintf('psi1 = %f, psi2 = %f, F = %f, CG iterations = %d, Avg line search iterations = %f\n', ...
            psi1_opt, psi2_opt, F_opt, cg_iter_count, avg_secant_iter);
        
        % 设置绘图范围
        psi1_vals = linspace(min(iter_points(:, 1))-1, max(iter_points(:, 1))+1, 50);
        psi2_vals = linspace(min(iter_points(:, 2))-1, max(iter_points(:, 2))+1, 50);

        % 计算 F 在 (psi1, psi2) 网格上的值
        [Psi1, Psi2] = meshgrid(psi1_vals, psi2_vals);
        F_vals = F(Psi1, Psi2);

        % 绘制自由能的二维颜色图
        figure;
        pcolor(Psi1, Psi2, F_vals);
        shading interp;  % 颜色平滑
        colorbar;        % 显示颜色条
        colormap turbo;
        hold on;

        % 添加等高线
        contour(Psi1, Psi2, F_vals, 20, 'LineColor', 'k', 'LineWidth', 0.5); % 20 条等高线

        % 添加标签
        xlabel('\psi_{CDW}');
        ylabel('\psi_{lattice}');
        title('free energy');

        % 绘制共轭梯度法的迭代轨迹
        plot(iter_points(:,1), iter_points(:,2), 'r-o', 'LineWidth', 0.5, 'MarkerSize', 4);
    end
end

% function [alpha, secant_iter] = line_search(F, grad_F, psi1, psi2, d)
%     % 使用修正的弦截法进行线搜索，确保方向导数的斜率为正并找到极小值点
%     % 同时记录每次调用时的弦截法迭代次数
% 
%     % 初始化步长值
%     alpha1 = 0;      % 初始步长
%     alpha2 = 0.1;    % 初始猜测步长
%     tol = 1e-8;      % 容差
%     max_iter = 100;   % 最大迭代次数
%     growth_factor = 2.0;  % 每次增加步长的倍数
%     alpha_max = 1;     % 步长的上限
%     secant_iter = 0;      % 记录弦截法的迭代次数
% 
%     % 计算初始方向导数
%     g1 = direction_derivative(F, grad_F, psi1, psi2, d, alpha1); % g(alpha1)
%     g2 = direction_derivative(F, grad_F, psi1, psi2, d, alpha2); % g(alpha2)
% 
%     for i = 1:max_iter
%         secant_iter = secant_iter + 1; % 记录弦截法的迭代次数
% 
%         % 计算方向导数的斜率
%         slope = (g2 - g1) / (alpha2 - alpha1);
%         
%         % 如果斜率为负，增加 alpha2 直到斜率为正
%         while slope < 0
%             alpha2 = alpha2 * growth_factor;  % 增加步长
%             g2 = direction_derivative(F, grad_F, psi1, psi2, d, alpha2); % 重新计算方向导数
%             slope = (g2 - g1) / (alpha2 - alpha1);  % 更新斜率
%         end
% 
%         % 使用弦截法更新步长
%         alpha_new = alpha2 - g2 * (alpha2 - alpha1) / (g2 - g1);
%         
%         % 限制步长为 alpha_max 以内
%         alpha_new = min(alpha_new, alpha_max);
% 
%         % 计算新的方向导数
%         g_new = direction_derivative(F, grad_F, psi1, psi2, d, alpha_new);
% 
%         % 检查收敛条件：步长的变化足够小，并且方向导数接近零
%         if abs(alpha_new - alpha2) < tol && abs(g_new) < tol
%             alpha = alpha_new;
%             return;
%         end
% 
%         % 更新步长和方向导数
%         alpha1 = alpha2;
%         g1 = g2;
%         alpha2 = alpha_new;
%         g2 = g_new;
%     end
%     
%     % 如果没有在最大迭代次数内收敛，返回最后的 alpha2
%     alpha = alpha2;
% end

function [alpha, bisect_iter] = line_search(F, grad_F, psi1, psi2, d)
    % 使用二分法进行线搜索，找到方向导数为零的步长 alpha
    % 确认初始方向导数为负，否则报错
    % 按恒定步长增加 alpha，直到方向导数变为正
    % 在变号区间内使用二分法求解零点

    % 初始化步长值
    alpha1 = 0;      % 初始步长
    alpha2 = 0.1;    % 初始猜测步长
    tol = 1e-6;      % 容差
    max_iter = 100;  % 最大迭代次数
    step_size = 0.1; % 恒定增加步长的大小
    bisect_iter = 0; % 记录二分法的迭代次数

    % 计算初始方向导数
    g1 = direction_derivative(F, grad_F, psi1, psi2, d, alpha1); % g(alpha1)
    
    % 确认初始方向导数是负的
    if g1 >= 0
        error('Initial direction derivative must be negative.');
    end

    % 增加 alpha，直到方向导数变为正
    g2 = direction_derivative(F, grad_F, psi1, psi2, d, alpha2); % g(alpha2)
    while g2 <= 0
        alpha1 = alpha2; % 更新下界
        alpha2 = alpha2 + step_size; % 增加 alpha
        g1 = g2;
        g2 = direction_derivative(F, grad_F, psi1, psi2, d, alpha2); % 重新计算方向导数
        
%         % 限制 alpha 的最大值
%         if alpha2 > 50  % 可以根据具体问题调整最大步长限制
%             error('Alpha too large, line search failed to find positive direction derivative.');
%         end
    end

    % 方向导数变号，进入二分法
    for i = 1:max_iter
        bisect_iter = bisect_iter + 1;

        % 中点
        alpha_mid = (alpha1 + alpha2) / 2;
        g_mid = direction_derivative(F, grad_F, psi1, psi2, d, alpha_mid);

        % 检查收敛条件
        if abs(g_mid) < tol || abs(alpha2 - alpha1) < tol
            alpha = alpha_mid;
            return;
        end

        % 更新区间：选择包含零点的区间
        if g_mid < 0
            alpha1 = alpha_mid;
        else
            alpha2 = alpha_mid;
        end
    end

    % 如果没有在最大迭代次数内收敛，返回中间值
    alpha = (alpha1 + alpha2) / 2;
end

% 计算给定步长 alpha 下的方向导数
function g = direction_derivative(F, grad_F, psi1, psi2, d, alpha)
    % psi 在步长 alpha 下的值
    psi1_new = psi1 + alpha * d(1);
    psi2_new = psi2 + alpha * d(2);
    
    % 计算在新点处的梯度
    grad_new = grad_F(psi1_new, psi2_new);
    
    % 计算方向导数 g(alpha) = grad(F) · d
    g = grad_new' * d;
end
