function out = plot_sepPsi_psiMax_bumpPsi_vs_mu_separate()
% Plot sep_psi, psi_max, bump_psi vs mu (your first column B) in 3 separate figures.

txt = [
"5.4  0        0        0        0"
"5.3  1.65062  3.78881  3.78881  0.469127"
"5.2  1.64752  3.78568  3.78568  0.494106"
"5.1  1.64752  3.97886  3.97886  0.519258"
"5.0  1.64924  4.30822  4.30822  0.547151"
"4.9  1.64924  4.42204  4.42204  0.575592"
"4.8  1.64736  4.06841  4.32631  0.607394"
"4.6  1.64874  4.59345  4.59345  0.670943"
"4.4  1.64874  4.76562  4.76562  0.730198"
"4.2  1.64503  3.99988  4.87830  0.772750"
"4.0  1.64628  4.18319  5.04432  0.794037"
"3.8  1.64753  4.27760  5.18830  0.798820"
"3.6  1.64884  4.46352  5.49519  0.777873"
"3.4  1.65027  4.85199  5.78344  0.765600"
"3.2  1.65188  5.34699  6.04031  0.788315"
"3.0  1.64865  3.79645  6.24846  0.816057"
"2.8  1.65078  4.83423  6.41687  0.829607"
"2.6  1.65323  5.68342  6.58111  0.825189"
"2.4  1.65321  5.58660  6.69819  0.77909"
"2.2  1.65090  4.39072  6.81214  0.712563"
"2.0  1.65229  4.87771  7.23034  0.801517"
"1.8  1.65080  3.89053  7.76062  0.92117"
"1.6  1.65470  5.20994  8.21721  1.01594"
"1.4  1.65352  4.55631  8.60090  1.11167"
"1.2  1.65237  3.88241  8.95898  1.15425"
"1.0  1.65648  4.80817  9.22145  1.14202"
"0.8  1.65523  4.23533  9.41635  1.12040"
"0.6  1.65384  3.68170  9.60868  1.09194"
"0.4  1.65763  4.39552  9.76030  1.13965"
"0.2  1.65583  3.88449  9.84675  1.21125"
"0.0  1.65970  4.58287  9.88432  1.21661"
];

M = parse_numeric_table_(txt);
mu = M(:,1);          % your x-axis (call it mu if you want)
sep_psi  = M(:,3);
psi_max  = M(:,4);
bump_psi = M(:,5);

% sort ascending for nicer lines
[mu, idx] = sort(mu, "ascend");
sep_psi  = sep_psi(idx);
psi_max  = psi_max(idx);
bump_psi = bump_psi(idx);

% ---- Figure 1: sep_psi ----
figure("Color","w","Units","pixels","Position",[120 120 820 480], "Name","sep_psi vs mu");
ax = axes();
plot(ax, mu, sep_psi, "o-","LineWidth",1.8,"MarkerSize",6);
grid(ax,"on"); box(ax,"on");
xlabel(ax,"mu (or polar\_mu)"); ylabel(ax,"sep\_psi");
title(ax,"sep\_psi vs mu","FontWeight","normal");

% ---- Figure 2: psi_max ----
figure("Color","w","Units","pixels","Position",[160 160 820 480], "Name","psi_max vs mu");
ax = axes();
plot(ax, mu, psi_max, "o-","LineWidth",1.8,"MarkerSize",6);
grid(ax,"on"); box(ax,"on");
xlabel(ax,"mu (or polar\_mu)"); ylabel(ax,"psi\_max");
title(ax,"psi\_max vs mu","FontWeight","normal");

% ---- Figure 3: bump_psi ----
figure("Color","w","Units","pixels","Position",[200 200 820 480], "Name","bump_psi vs mu");
ax = axes();
plot(ax, mu, bump_psi, "o-","LineWidth",1.8,"MarkerSize",6);
grid(ax,"on"); box(ax,"on");
xlabel(ax,"mu (or polar\_mu)"); ylabel(ax,"bump\_psi");
title(ax,"bump\_psi vs mu","FontWeight","normal");

out = struct();
out.mu = mu;
out.sep_psi = sep_psi;
out.psi_max = psi_max;
out.bump_psi = bump_psi;
end

function M = parse_numeric_table_(lines)
n = numel(lines);
M = nan(n,5);
for i = 1:n
    s = strtrim(lines(i));
    s = regexprep(s, "\s+", " ");
    v = sscanf(s, "%f");
    if numel(v) < 5
        error("Line %d has <5 numbers: %s", i, s);
    end
    M(i,:) = v(1:5).';
end
end