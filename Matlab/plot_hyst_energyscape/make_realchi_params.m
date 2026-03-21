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

    % =================== q-point selection ===================
    par.iq_pick = -267;
    par.jq_pick = -267;

    % =================== chi -> coefficients =================
    par.invV       = -10.14;
    par.chi_scaler = -500;
    par.T_target   = 6.5;          % used only for picking one T from loaded list
    par.T_tol      = 1e-8;

    % =================== coupling ============================
    par.lambda0 = 2;
    par.lambda  = 2;

    % =================== psi window ==========================
    par.psi1_lim = [-12, 12];
    par.psi2_lim = [-8, 8];

    % =================== doping reference ====================
    par.dop_span = 1.0;

    % =================== define pt1(dop1, T1), pt2(dop2, T2) on the lattice trans boundary ==============
    % =================== T = @(dop) -A * dop^2 + C for lattice mode boundary (A,C detemined in coeff.m) ====================
    par.Lat_pt1 = [-3, 5];
    par.Lat_pt2 = [-2.7, 8];
    
    
    % =================== a2 = alpha * (T - Tc) is predefined.(alpha > 0) ====================
    par.Lat_Tc = 10;
    par.Lat_alpha = 5e-3;
    

    % =================== artificial magnetic field (meV) =========
    par.artificial_polar = 0; 

    % =========================================================
    % =================== IO + FILTERS ========================
    % =========================================================
    par.io = struct();
    par.io.default_root = "D:\OneDrive - Emory\Rhombohedral_SC\rhombohedral_project\data\";

    % -------- NEW: range filters --------
    % Temperature range is parsed from folder name (preferred).
    % Set [-Inf Inf] to disable.
    par.io.T_range_K = [0, 100];
    par.io.mu_range  = [0, 100];
    par.io.dop_range = [0, 100];

    % =========================================================
    % =================== HYSTERESIS SWEEP ====================
    % =========================================================
    par.hyst = struct();

    % scanning direction
    par.hyst.forward_direction = "ascend";   % "ascend" or "descend"

    % --- grid type: choose ONE ---
    par.hyst.grid_mode = "N";        % "N" or "step"
    par.hyst.N         = 200;        % used if grid_mode="N"
    par.hyst.step      = 0.002;      % used if grid_mode="step" (unit = your doping unit)

    % initial seed for continuation
    par.hyst.psi1_0 = 1e-4;
    par.hyst.psi2_0 = 0;

    % (optional) a marker doping (not required)
    par.hyst.turn_dop = 1.02;

    % =========================================================
    % ============ Dense minima search (robust) ================
    % =========================================================
    par.min = struct();

    par.min.Npsi1 = 81;
    par.min.Npsi2 = 81;

    par.min.use_gaussian_smooth = true;
    par.min.smooth_sigma        = 0.7;

    par.min.keep_topK = 40;

    par.min.seed_jitter  = 0.12;
    par.min.seed_repeats = 2;
    par.min.max_iter     = 450;
    par.min.tol_fun      = 1e-10;
    par.min.tol_x        = 1e-8;

    par.min.cluster_tol  = 3e-2;
    par.min.fd_h         = 2e-3;
    par.min.min_eig_eps  = 1e-6;

    par.min.force_origin_and_nearby = true;
    par.min.nearby_delta = 0.2;

    % =========================================================
    % =================== PLOTTING SMOOTH =====================
    % =========================================================
    par.smooth = struct();

    % Lorentz broadening for psi(doping) plots ONLY
    par.smooth.use_lorentz = false;

    % gamma in SAME UNIT as doping axis.
    % Here we keep your convention: gamma_factor * median(diff(dop_grid)).
    %   gamma_factor <0 or NaN -> OFF
    %   gamma_factor = 0       -> auto (1.5 * dx)
    par.smooth.lorentz_gamma_factor = 2;  % was "lorentz_gamma_dop" in your old code

    % optional: also smooth chi(dop) in panel 1
    par.smooth.smooth_chi_curve = false;

    % =========================================================
    % =================== plotting style ======================
    % =========================================================
    par.plot = struct();
    par.plot.fontSize = 14;
    par.plot.show_gamma_in_title = false;
    par.plot.save_dir_name = "plot";
    par.plot.export_dpi = 300;


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