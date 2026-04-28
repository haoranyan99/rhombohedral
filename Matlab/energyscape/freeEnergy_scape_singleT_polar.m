function S = build_shifted_state_singleT_(R, U_pick, u_tag, polar_meV)

    S = struct('ok',false,'msg',"",'chi_used',NaN,'dop_eval',NaN, ...
               'mu_used',NaN,'dmu_eV',NaN);

    % -------------------------------------------------
    % Use ALL raw points in this selected T folder
    % -------------------------------------------------
    mu0  = R.muVals(:);
    dop0 = R.dopVals(:);
    chi0 = R.chiVals(:);

    ok = isfinite(mu0) & isfinite(dop0) & isfinite(chi0);
    mu0  = mu0(ok);
    dop0 = dop0(ok);
    chi0 = chi0(ok);

    if numel(mu0) < 3
        S.msg = 'Not enough raw chi points in selected T folder.';
        return;
    end

    % collapse duplicate mu over the whole T folder
    [mu0, chi0, dop0] = collapse_same_mu_(mu0, chi0, dop0);

    if numel(mu0) < 3
        S.msg = 'Too few unique mu points after collapsing duplicates.';
        return;
    end

    % -------------------------------------------------
    % choose center mu from current UI selection
    % -------------------------------------------------
    if strcmpi(u_tag, 'mu')
        mu_center = double(U_pick);
    else
        % current slider is folder doping -> choose the mu whose doping is closest
        [~, id0] = min(abs(dop0 - U_pick));
        mu_center = mu0(id0);
    end

    dmu = 0.5 * double(polar_meV) * 1e-3;
    S.dmu_eV = dmu;

    % -------------------------------------------------
    % evaluate shifted chi using nearest left/right points
    % -------------------------------------------------
    if abs(dmu) < 1e-18
        chi_used = interp_or_neighbor_mean_(mu0, chi0, mu_center);
    else
        chi_p = interp_or_neighbor_mean_(mu0, chi0, mu_center + dmu);
        chi_m = interp_or_neighbor_mean_(mu0, chi0, mu_center - dmu);
        chi_used = 0.5 * (chi_p + chi_m);
    end

    % doping evaluated at center mu
    dop_eval = interp_or_neighbor_mean_(mu0, dop0, mu_center);

    if ~(isfinite(chi_used) && isfinite(dop_eval) && isfinite(mu_center))
        S.msg = 'Shift-averaged chi is not finite.';
        return;
    end

    S.ok = true;
    S.chi_used = chi_used;
    S.dop_eval = dop_eval;
    S.mu_used  = mu_center;
end

function y = interp_or_neighbor_mean_(x, v, xq)
% For query xq:
% 1) if xq lies inside data range, use linear interpolation
% 2) if interpolation is unavailable / outside range:
%       use mean of nearest left and right points
% 3) if only one side exists, use that nearest point

    x = double(x(:));
    v = double(v(:));

    ok = isfinite(x) & isfinite(v);
    x = x(ok);
    v = v(ok);

    if isempty(x)
        y = NaN;
        return;
    end

    [x, ord] = sort(x);
    v = v(ord);

    % unique x only
    [xu, ia] = unique(x);
    vu = v(ia);

    if numel(xu) == 1
        y = vu(1);
        return;
    end

    % first try linear interpolation inside range
    if xq >= xu(1) && xq <= xu(end)
        y = interp1(xu, vu, xq, 'linear');
        if isfinite(y)
            return;
        end
    end

    % otherwise use nearest left/right mean
    iL = find(xu <= xq, 1, 'last');
    iR = find(xu >= xq, 1, 'first');

    if ~isempty(iL) && ~isempty(iR)
        if iL == iR
            y = vu(iL);
        else
            y = 0.5 * (vu(iL) + vu(iR));
        end
    elseif ~isempty(iL)
        y = vu(iL);
    elseif ~isempty(iR)
        y = vu(iR);
    else
        y = NaN;
    end
end

function [mu2, chi2, dop2] = collapse_same_mu_(mu, chi, dop)
    mu  = double(mu(:));
    chi = double(chi(:));
    dop = double(dop(:));

    ok = isfinite(mu) & isfinite(chi) & isfinite(dop);
    mu  = mu(ok);
    chi = chi(ok);
    dop = dop(ok);

    [mu_u,~,ic] = unique(mu);
    chi_m = accumarray(ic, chi, [], @mean);
    dop_m = accumarray(ic, dop, [], @mean);

    [mu2, ord] = sort(mu_u);
    chi2 = chi_m(ord);
    dop2 = dop_m(ord);
end