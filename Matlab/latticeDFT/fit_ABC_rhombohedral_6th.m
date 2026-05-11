function fit_ABC_rhombohedral_6th()
clc; close all;

% =========================
% Digitized data
% =========================
x = [ ...
    0.00
    0.10
    0.20
    0.30
    0.40
    0.52
    0.63
    0.74
    0.85
%     0.96
%     1.05
%     1.15
%     1.26
%     1.37
%     1.48
%     1.59
%     1.69
%     1.79
%     1.90
%     2.00 
    ];

y_ABC = [ ...
      0
    180
    520
    675
    560
    435
    510
    700
    785
%     695
%     615
%     690
%     820
%     850
%     750
%     660
%     690
%     750
%     715
%     580 
    ];

% =========================
% 6th-order polynomial fit
% =========================
poly_order = 6;
p = polyfit(x, y_ABC, poly_order);

x_fit = linspace(min(x), max(x), 800);
y_fit = polyval(p, x_fit);

% fitted value at original points
y_pred = polyval(p, x);

% R^2
SS_res = sum((y_ABC - y_pred).^2);
SS_tot = sum((y_ABC - mean(y_ABC)).^2);
R2 = 1 - SS_res / SS_tot;

% =========================
% Print fit function
% =========================
fprintf('\n6th-order polynomial fit:\n');
fprintf('E(x) = %.10g*x^6 + %.10g*x^5 + %.10g*x^4 + %.10g*x^3 + %.10g*x^2 + %.10g*x + %.10g\n', p);
fprintf('R^2 = %.6f\n\n', R2);

% =========================
% Save data
% =========================
fit_data.x = x;
fit_data.y_ABC = y_ABC;
fit_data.p6 = p;
fit_data.R2 = R2;
fit_data.x_fit = x_fit;
fit_data.y_fit = y_fit;

% =========================
% Plot
% =========================
figure('Color','w');
hold on;

plot(x, y_ABC, 'o', ...
    'LineWidth',2.2, ...
    'MarkerSize',9, ...
    'DisplayName','Digitized ABC data');

plot(x_fit, y_fit, '-', ...
    'LineWidth',3.0, ...
    'DisplayName','6th-order fit');

xlabel('$X/a$', 'Interpreter','latex', 'FontSize',24);
ylabel('$E_{\mathrm{inter}}\;(\mathrm{meV})$', ...
    'Interpreter','latex', 'FontSize',24);

legend('Interpreter','latex', ...
    'FontSize',18, ...
    'Location','best');

set(gca, ...
    'FontSize',20, ...
    'LineWidth',1.8, ...
    'TickDir','in', ...
    'Box','on');

xlim([0 2.05]);
ylim([0 900]);
grid on;

end