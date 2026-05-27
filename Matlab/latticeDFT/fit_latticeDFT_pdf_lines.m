function fit_latticeDFT_pdf_lines()
clc; close all;

% Data extracted from vector paths in 6L_rhombo_1v5_2v4_MK80.pdf.
% x unit: X/a. y unit: eV.
x_blue = [ ...
    0.0000000000
    0.1052845547
    0.2105691032
    0.3158536579
    0.4211382126
    0.5264227611
    0.6317073158
    0.7369918705
    0.8422764190
    0.9475609737
    1.0524390202
    1.1577235749
    1.2630081295
    1.3682926842
    1.4735772328
    1.5788617874
    1.6841463421
    1.7894308907
    1.8947154453
    2.0000000000 ];

y_blue_eV = [ ...
    0.0000000000
    8.2479379290
    23.8854213836
    30.9592335085
    25.6355020804
    19.5591469610
    22.9636947209
    31.9740112643
    35.3354341270
    31.0984808851
    27.7621535815
    31.5381605044
    37.9272058609
    38.5923122629
    33.4626671471
    29.7705710748
    31.5960488194
    34.6726592773
    32.3495077639
    25.6095291023 ] / 40.0;

x_red = [ ...
    0.0000000000
    0.1052845547
    0.2105691032
    0.3158536579
    0.6317073158
    0.7369918705
    0.9475609737
    1.0524390202
    1.1577235749
    1.2630081295
    1.3682926842
    1.5788617874
    1.6841463421
    1.7894308907
    1.8947154453
    2.0000000000 ];

y_red_eV = [ ...
    0.0000000000
    4.1192107769
    11.9216967350
    15.4016659016
    11.2663719667
    15.7570366987
    15.1572989316
    13.4210493005
    15.2202474827
    18.3638461803
    18.6093990014
    13.9783712620
    14.7839004763
    16.1643477404
    14.8591051740
    11.2586408357 ] / 40.0;

poly_order = 6;

% Fit only the X/a <= 1 region.
fit_idx_blue = find(x_blue <= 0.8);
fit_idx_red = find(x_red <= 0.8);

x_blue_fit_data = x_blue(fit_idx_blue);
y_blue_fit_data = y_blue_eV(fit_idx_blue);
x_red_fit_data = x_red(fit_idx_red);
y_red_fit_data = y_red_eV(fit_idx_red);

[p_blue, R2_blue] = local_polyfit_r2( ...
    x_blue_fit_data, y_blue_fit_data, poly_order);

[p_red, R2_red] = local_polyfit_r2( ...
    x_red_fit_data, y_red_fit_data, poly_order);

x_fit_blue = linspace(min(x_blue_fit_data), max(x_blue_fit_data), 600);
x_fit_red = linspace(min(x_red_fit_data), max(x_red_fit_data), 600);
y_blue_fit = polyval(p_blue, x_fit_blue);
y_red_fit = polyval(p_red, x_fit_red);

fprintf('\nBlue line %dth-order fit, y in eV:\n', poly_order);
local_print_poly('E_blue', p_blue);
fprintf('R^2_blue = %.8f\n', R2_blue);

fprintf('\nRed line %dth-order fit, y in eV:\n', poly_order);
local_print_poly('E_red', p_red);
fprintf('R^2_red = %.8f\n\n', R2_red);

fit_data = struct();
fit_data.source_pdf = '6L_rhombo_1v5_2v4_MK80.pdf';
fit_data.poly_order = poly_order;
fit_data.x_blue = x_blue;
fit_data.y_blue_eV = y_blue_eV;
fit_data.x_red = x_red;
fit_data.y_red_eV = y_red_eV;
fit_data.fit_idx_blue = fit_idx_blue;
fit_data.fit_idx_red = fit_idx_red;
fit_data.x_blue_fit_data = x_blue_fit_data;
fit_data.y_blue_fit_data_eV = y_blue_fit_data;
fit_data.x_red_fit_data = x_red_fit_data;
fit_data.y_red_fit_data_eV = y_red_fit_data;
fit_data.p_blue = p_blue;
fit_data.p_red = p_red;
fit_data.R2_blue = R2_blue;
fit_data.R2_red = R2_red;
fit_data.x_fit_blue = x_fit_blue;
fit_data.x_fit_red = x_fit_red;
fit_data.y_blue_fit_eV = y_blue_fit;
fit_data.y_red_fit_eV = y_red_fit;

save('latticeDFT_pdf_line_fit.mat', 'fit_data');

tbl_blue = table(x_blue, y_blue_eV, 'VariableNames', {'x_over_a', 'E_eV'});
tbl_red = table(x_red, y_red_eV, 'VariableNames', {'x_over_a', 'E_eV'});
writetable(tbl_blue, 'latticeDFT_pdf_blue_line_points.csv');
writetable(tbl_red, 'latticeDFT_pdf_red_line_points.csv');

figure('Color', 'w', 'Position', [200 200 390 280]);
hold on;

plot(x_blue, y_blue_eV, '-o', ...
    'Color', [0.1215686275 0.2274509804 0.4078431373], ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', 5.5, ...
    'LineWidth', 1.3, ...
    'DisplayName', 'DFT simulation');

plot(x_fit_blue, y_blue_fit, 'r--', ...
    'LineWidth', 1.4, ...
    'DisplayName', 'blue 6th-order potential');

plot(x_red, y_red_eV, '-s', ...
    'Color', [0.7529411765 0.2235294118 0.168627451], ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', 5.0, ...
    'LineWidth', 1.3, ...
    'DisplayName', 'red DFT simulation');

plot(x_fit_red, y_red_fit, '--', ...
    'Color', [0.65 0.1 0.1], ...
    'LineWidth', 1.4, ...
    'DisplayName', 'red 6th-order potential');

xlabel('lattice distortion (X/a)', 'FontSize', 12);
ylabel('energy [eV]', 'FontSize', 12);

legend('Location', 'southeast', 'FontSize', 10, 'Box', 'on');
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'in', 'Box', 'on');
xlim([0 2]);
ylim([0 1]);
text(0.03, 0.93, '(b)', 'FontSize', 12);

savefig(gcf, 'latticeDFT_pdf_line_fit.fig');

end

function [p, R2] = local_polyfit_r2(x, y, order)
p = polyfit(x(:), y(:), order);
y_pred = polyval(p, x(:));
SS_res = sum((y(:) - y_pred).^2);
SS_tot = sum((y(:) - mean(y(:))).^2);
R2 = 1 - SS_res / SS_tot;
end

function local_print_poly(name, p)
order = numel(p) - 1;
fprintf('%s(x) = ', name);
for i = 1:numel(p)
    pow = order - i + 1;
    coeff = p(i);
    if i == 1
        fprintf('%.12g*x^%d', coeff, pow);
    else
        fprintf(' %+.12g*x^%d', coeff, pow);
    end
end
fprintf('\n');
end
