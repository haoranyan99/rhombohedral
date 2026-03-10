function par = make_realchi_params(noAsk)
% make_realchi_params
% Centralized params for hysteresis (real chi) + plotting + smoothing.
%
% Usage:
%   par = make_realchi_params();         % default: ASK to backup
%   par = make_realchi_params(true);     % no UI, no backup prompt

    if nargin < 1
        noAsk = false;
    end

    par = struct();

    par.io.T_range_K   = [3.999, 4.001];     % 温度范围（K），从文件夹名解析
    par.io.dop_range   = [0.90, 1.10];       % doping范围（用header读出后过滤）

    % =================== q-point selection ===================
    par.iq_pick = -133;
    par.jq_pick = -133;

    % =================== chi -> coefficients =================
    par.invV       = -10.21;
    par.chi_scaler = -500;
    par.T_target   = 4;          % target temperature to pick nearest T in data
    par.T_tol      = 1e-8;

    % =================== coupling ============================
    par.lambda0 = 4;
    par.lambda  = 4;

    % =================== psi window ==========================
    par.psi1_lim = [-12, 12];
    par.psi2_lim = [-10, 10];

    % =================== doping reference ====================
    par.dop_span = 1.0;

    % =================== coefficient: b2(T) ==================
    par.Tc     = 20;
    par.Tslope = 2.5 / (par.Tc - par.T_target);

    % =================== coefficient: a2(T,dop) ==============
    par.a2_T = 4.044e-3;
    par.a2_D = 5;
    par.a2_c = 5.71;

    % =========================================================
    % =================== HYSTERESIS SWEEP ====================
    % =========================================================
    par.hyst = struct();

    % --- grid type: choose ONE ---
    par.hyst.grid_mode = "N";        % "N" or "step"
    par.hyst.N         = 200;        % used if grid_mode="N"
    par.hyst.step      = 0.002;      % used if grid_mode="step" (unit = your doping unit)

    % initial seed for continuation
    par.hyst.psi1_0 = 1e-4;
    par.hyst.psi2_0 = 0;

    % (optional) a marker doping (you used before; not necessary)
    par.hyst.turn_dop = 1.02;

    % =========================================================
    % ============ Dense minima search (robust) ================
    % =========================================================
    par.min = struct();

    % coarse grid for dense scan
    par.min.Npsi1 = 81;
    par.min.Npsi2 = 81;

    % smooth Fgrid for seed finding (this is your old behavior)
    par.min.use_gaussian_smooth = true;
    par.min.smooth_sigma        = 0.7;

    % keep best seeds
    par.min.keep_topK = 40;

    % fminsearch details
    par.min.seed_jitter  = 0.12;
    par.min.seed_repeats = 2;
    par.min.max_iter     = 450;
    par.min.tol_fun      = 1e-10;
    par.min.tol_x        = 1e-8;

    % clustering + Hessian filter
    par.min.cluster_tol  = 3e-2;
    par.min.fd_h         = 2e-3;
    par.min.min_eig_eps  = 1e-6;

    % add origin + nearby seeds
    par.min.force_origin_and_nearby = true;
    par.min.nearby_delta = 0.2;

    % =========================================================
    % =================== PLOTTING SMOOTH =====================
    % =========================================================
    par.smooth = struct();

    % Lorentz broadening for psi(doping) plots ONLY
    par.smooth.use_lorentz = true;

    % gamma in SAME UNIT as doping axis.
    % recommendation: gamma ~ (1~2)*median(diff(dop_grid))
    par.smooth.lorentz_gamma_dop = 3;  % 0 -> auto; negative/NaN -> off

    % choose whether also smooth chi(dop) curve shown in panel 1
    par.smooth.smooth_chi_curve = false;

    % =========================================================
    % =================== plotting style ======================
    % =========================================================
    par.plot = struct();
    par.plot.fontSize = 14;
    par.plot.show_gamma_in_title = false;
    par.plot.save_dir_name = "plot";
    par.plot.export_dpi = 300;

    % default root for selecting folder
    par.io = struct();
    par.io.default_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\";

    % =========================================================
    % =================== optional backup =====================
    % =========================================================
    if noAsk
        return;
    end

    choice = questdlg('Save parameter backup?', ...
                      'Parameter Backup', ...
                      'Yes','No','No');

    if strcmp(choice,'Yes')
        save_dir = uigetdir(pwd, 'Select folder to save parameter backup');
        if save_dir ~= 0
            ts = datestr(now,'yyyy_mm_dd_HHMMSS');
            matfile = fullfile(save_dir, sprintf('realchi_params_%s.mat', ts));
            save(matfile, 'par');

            jsonfile = fullfile(save_dir, sprintf('realchi_params_%s.json', ts));
            json_text = jsonencode(par, 'PrettyPrint', true);
            fid = fopen(jsonfile, 'w');
            fwrite(fid, json_text, 'char');
            fclose(fid);

            fprintf('[param backup] saved:\n  %s\n  %s\n', matfile, jsonfile);
        end
    end
end