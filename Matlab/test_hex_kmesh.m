function test_hex_kmesh_RG()

clc; close all;

% =========================
% Parameters
% =========================
n_k = 5;
dk_frac = 1.0;

a = 2.46;

% =========================
% Real-space primitive vectors
% Same as RG_Structure.h
% a1 = (0, -a)
% a2 = (sqrt(3)/2*a, 1/2*a)
% =========================
a1 = [0.0, -a];
a2 = [sqrt(3)/2 * a, 0.5 * a];

% =========================
% Reciprocal vectors
% Same as C++:
% b1 = 2pi * ( a2_y, -a2_x) / area
% b2 = 2pi * (-a1_y,  a1_x) / area
% =========================
area = a1(1)*a2(2) - a1(2)*a2(1);
twoPi = 2*pi;

b1 = twoPi * [ a2(2), -a2(1)] / area;
b2 = twoPi * [-a1(2),  a1(1)] / area;

Gamma = [0, 0];
M  = 0.5 * b1;
K  = (b1 + b2) / 3;
Kp = 2 * (b1 + b2) / 3;

dk1 = dk_frac * b1;
dk2 = dk_frac * b2;

% =========================
% Generate hex mesh
% Same as new C++ logic
% =========================
kx = [];
ky = [];
iq_list = [];
jq_list = [];

for iq = -n_k:n_k

    jq_min = max(-n_k, -iq - n_k);
    jq_max = min( n_k, -iq + n_k);

    fprintf('iq = %2d : jq = [%2d, %2d]\n', ...
        iq, jq_min, jq_max);

    for jq = jq_min:jq_max

        k = K + iq * dk1 + jq * dk2;

        kx(end+1) = k(1);
        ky(end+1) = k(2);

        iq_list(end+1) = iq;
        jq_list(end+1) = jq;
    end
end

% =========================
% Plot
% =========================
figure;
hold on; box on;

scatter(kx, ky, 80, 'filled');

for i = 1:numel(kx)
    text(kx(i), ky(i), ...
        sprintf('(%d,%d)', iq_list(i), jq_list(i)), ...
        'FontSize', 8, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
end

% plot high-symmetry points
scatter(Gamma(1), Gamma(2), 120, 'd', 'filled');
text(Gamma(1), Gamma(2), '\Gamma', 'FontSize', 14);

scatter(M(1), M(2), 120, 's', 'filled');
text(M(1), M(2), 'M', 'FontSize', 14);

scatter(K(1), K(2), 120, '^', 'filled');
text(K(1), K(2), 'K', 'FontSize', 14);

scatter(Kp(1), Kp(2), 120, 'v', 'filled');
text(Kp(1), Kp(2), 'K''', 'FontSize', 14);

% plot reciprocal basis vectors from K
quiver(K(1), K(2), dk1(1), dk1(2), 0, 'LineWidth', 2);
text(K(1)+dk1(1), K(2)+dk1(2), 'dk_1 = dk\_frac b_1');

quiver(K(1), K(2), dk2(1), dk2(2), 0, 'LineWidth', 2);
text(K(1)+dk2(1), K(2)+dk2(2), 'dk_2 = dk\_frac b_2');

axis equal;
grid on;

xlabel('k_x');
ylabel('k_y');

title(sprintf('RG local K hex mesh, n_k = %d, dk\\_frac = %.4g', ...
    n_k, dk_frac));

fprintf('\nExpected points = %d\n', 1 + 3*n_k*(n_k+1));
fprintf('Actual points   = %d\n', numel(kx));

fprintf('\na1 = (%.8f, %.8f)\n', a1(1), a1(2));
fprintf('a2 = (%.8f, %.8f)\n', a2(1), a2(2));
fprintf('b1 = (%.8f, %.8f)\n', b1(1), b1(2));
fprintf('b2 = (%.8f, %.8f)\n', b2(1), b2(2));
fprintf('K  = (%.8f, %.8f)\n', K(1), K(2));

end