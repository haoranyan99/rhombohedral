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
    par.iq_pick = 3;
    par.jq_pick = -3;

    % =================== chi -> coefficients =================
    par.invV       = 12.28;
    par.chi_scaler = 300; 

    par.T_target   = 2;          % used only for picking one T from loaded list
    par.T_tol      = 1e-8;

    % =================== coupling ============================
    par.lambda0 = 2;
    par.lambda  = 2;

    % =================== psi window ==========================
    par.psi1_lim = [-12, 12];
    par.psi2_lim = [-10, 10];

    % =================== doping reference ====================
    par.dop_span = 1.0;

    % =================== lattice transition boundary =========
    % The lattice X double-well boundary is defined by
    %   delta = (5/6)*b2^2 - a2*c2 = 0.
    %
    % "parabola" keeps the old T_boundary(dop) = -A*dop^2 + C.
    % "sigmoid_step" defines a monotonic boundary x_c(T):
    %   T = 0..T_knee_K  : slowly moves right
    %   T > T_knee_K    : rapidly approaches dop_right
    % The coefficient code inverts x_c(T) numerically to get T_boundary(dop).
    par.Lat_boundary_mode = "sigmoid_step";

    % Old parabola parameters, kept for quick comparison.
    par.Lat_pt1 = [-2, 2];
    par.Lat_pt2 = [0, 10];

    par.Lat_step = struct();
    par.Lat_step.T_min_K = 0;
    par.Lat_step.T_knee_K = 8;
    par.Lat_step.T_max_K = 10;
    par.Lat_step.dop_left = -2.0;
    par.Lat_step.dop_knee = -0.60;
    par.Lat_step.dop_right = 0.0;
    par.Lat_step.slow_power = 1.45;
    par.Lat_step.fast_width_K = 0.35;
    par.Lat_step.n_grid = 5001;

    
    % =================== a2 = alpha * (T - Tc) is predefined.(alpha > 0) ====================
    par.Lat_Tc = 0;
    par.Lat_alpha = 0.2;
    

    % =================== artificial magnetic field (meV) =========
    par.artificial_polar = 0; 

    % =========================================================
    % =================== IO + FILTERS ========================
    % =========================================================
    par.io = struct();
    par.io.default_root = "E:/rg_master/data";

    % -------- NEW: range filters --------
    % Temperature range is parsed from folder name (preferred).
    % Set [-Inf Inf] to disable.
    par.io.T_range_K = [0, 100];
    par.io.mu_range  = [0, 100];
    par.io.dop_range = [-Inf Inf];

    % =========================================================
    % =================== HYSTERESIS SWEEP ====================
    % =========================================================
    par.plot_paths = ["A","B","C","D","E"];
    par.plot_paths = ["A","B"];
    par.hyst = struct();

    % scanning direction
    par.hyst.forward_direction = "ascend";   % "ascend" or "descend"
    par.hyst.doping_bump = -1.44;
    

    % --- NEW: scan window on doping axis ---------------------
    % If NaN, use full available doping range from data.
    % If given, actual scan range will be:
    %   [max(min(dop0), min(scan_start,scan_end)), ...
    %    min(max(dop0), max(scan_start,scan_end))]
    par.hyst.scan_start = -2.5;
    par.hyst.scan_end   = -0.175;

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
