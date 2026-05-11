% =========================
% Digitized data from figure
% =========================

x_over_a = [ ...
    0.00
    0.10
    0.20
    0.30
    0.40
    0.52
    0.63
    0.74
    0.85
    0.96
    1.05
    1.15
    1.26
    1.37
    1.48
    1.59
    1.69
    1.79
    1.90
    2.00 ];

ABC_rhombohedral_meV = [ ...
      0
    180
    520
    675
    560
    435
    510
    700
    785
    695
    615
    690
    820
    850
    750
    660
    690
    750
    715
    580 ];

AB_Bernal_meV = [ ...
      0
    125
    390
    500
    400
    310
    390
    530
    575
    500
    455
    525
    625
    620
    520
    460
    495
    530
    460
    325 ];

% =========================
% Save MAT
% =========================
data.x_over_a = x_over_a;
data.ABC_rhombohedral_meV = ABC_rhombohedral_meV;
data.AB_Bernal_meV = AB_Bernal_meV;

save('LA_pulse_N40_E_inter_digitized.mat', ...
    'x_over_a', ...
    'ABC_rhombohedral_meV', ...
    'AB_Bernal_meV', ...
    'data');

% =========================
% Plot reproduction
% =========================
figure('Color','w');
hold on;

plot(x_over_a, ABC_rhombohedral_meV, ...
    '-o', ...
    'LineWidth',2.8, ...
    'MarkerSize',10);

plot(x_over_a, AB_Bernal_meV, ...
    '-s', ...
    'LineWidth',2.8, ...
    'MarkerSize',10);

xlabel('$X/a$', ...
    'Interpreter','latex', ...
    'FontSize',24);

ylabel('$E_{\mathrm{inter}}=\Delta E_{6\mathrm{L}}-2\Delta E_{3\mathrm{L}}\ (\mathrm{meV})$', ...
    'Interpreter','latex', ...
    'FontSize',24);

legend({'ABC rhombohedral','AB Bernal'}, ...
    'Interpreter','latex', ...
    'FontSize',20, ...
    'Location','northeast');

set(gca, ...
    'FontSize',20, ...
    'LineWidth',1.8, ...
    'TickDir','in', ...
    'Box','on');

xlim([0 2.05]);

grid on;