function plot_lattice_phase_boundary()
% Plot the lattice X-mode phase boundary implied by make_realchi_params.

this_dir = fileparts(mfilename("fullpath"));
addpath(this_dir, "-begin");

par = make_realchi_params(true);
coef = make_realchi_coeff(par);

dop = linspace(-2.2, 0.1, 500);
T = linspace(0, 10, 400);
[DOP, TT] = meshgrid(dop, T);

a2 = coef.a2(TT);
b2 = coef.b2(DOP);
c2 = coef.c2;
Delta = (5/6) .* b2.^2 - a2 .* c2;
Tbd = coef.T_boundary(dop);

figure("Name", "Lattice phase diagram", "Color", "w");
imagesc(dop, T, Delta < 0);
set(gca, "YDir", "normal");
hold on;
plot(dop, Tbd, "k-", "LineWidth", 2);

colormap([0.90 0.18 0.14; 0.22 0.45 0.88]);
cb = colorbar;
cb.Ticks = [0.25, 0.75];
cb.TickLabels = {'\Delta > 0', '\Delta < 0'};

xlabel("doping");
ylabel("T");
title("Lattice phase boundary");
xlim([min(dop), max(dop)]);
ylim([min(T), max(T)]);
box on;

end
