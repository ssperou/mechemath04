function plottingDP_GlobalTruncation()
    %% Dormand-Prince tableau 
    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;
                       5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
                       1/5, 0, 0, 0, 0, 0, 0;
                       3/40, 9/40, 0, 0, 0, 0, 0;
                       44/45, -56/15, 32/9, 0, 0, 0, 0;
                       19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0;
                       9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, 0;
                       35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
    DP_fixed = DormandPrince;
    DP_fixed.B = DormandPrince.B(1, :);

    %% Problem setup
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    V0 = [8; 0; 0; 1.5];
    t_span = [0 2*pi];
    X0 = V0(:);
    rate_func = @(t, V) gravity_rate_func(t, V, orbit_params);

    x_real = compute_planetary_motion(t_span(end), V0, orbit_params);
    x_real = x_real(:);

    %% Sweep settings
    n_pts = 40;
    p_DP = 5;
    growth_limit = 8;
    h_ref = 1e-2;                          % starting step guess for the adaptive method

    % Adaptive: sweep the tolerance
    e_list = logspace(-8, -2, n_pts);
    e_list = e_list(end:-1:1);

    % Fixed: sweep the step size
    h_list_fixed = logspace(-0.3, -2.5, n_pts);

    % Fit only where the error is in the clean power-law region:
    % above the roundoff floor, below the huge-step saturation.
    err_bounds_fixed    = [1e-10, 1e-3];
    err_bounds_adaptive = [1e-10, 1e-3];

    %% Fixed-step sweep
    error_fix = NaN(n_pts, 1);
    h_avg_fix = NaN(n_pts, 1);
    evals_fix = NaN(n_pts, 1);
    for i = 1:n_pts
        try
            [~, X_list, h_avg, evals] = explicit_RK_fixed_step_integration( ...
                rate_func, t_span, X0, h_list_fixed(i), DP_fixed);
            error_fix(i) = norm(X_list(end, :).' - x_real);
            h_avg_fix(i) = h_avg;
            evals_fix(i) = evals;
        catch ME
            warning('Fixed step failed at h = %g: %s', h_list_fixed(i), ME.message);
        end
    end

    %% Adaptive sweep
    error_var = NaN(n_pts, 1);
    h_avg_var = NaN(n_pts, 1);
    evals_var = NaN(n_pts, 1);
    fail_var  = NaN(n_pts, 1);
    for i = 1:n_pts
        try
            [~, X_list, h_avg, evals, failed] = explicit_RK_variable_step_integration( ...
                rate_func, t_span, X0, h_ref, DormandPrince, p_DP, e_list(i), growth_limit);
            error_var(i) = norm(X_list(end, :).' - x_real);
            h_avg_var(i) = h_avg;
            evals_var(i) = evals;
            fail_var(i)  = failed;
        catch ME
            warning('Adaptive step failed at tol = %g: %s', e_list(i), ME.message);
        end
    end

    %% Fits
    mask_fix = isfinite(error_fix) & isfinite(h_avg_fix) & ...
               error_fix >= err_bounds_fixed(1) & error_fix <= err_bounds_fixed(2);
    mask_var = isfinite(error_var) & isfinite(h_avg_var) & ...
               error_var >= err_bounds_adaptive(1) & error_var <= err_bounds_adaptive(2);

    [slope_fix, yfit_fix] = fit_loglog(h_avg_fix, error_fix, mask_fix);
    [slope_var, yfit_var] = fit_loglog(h_avg_var, error_var, mask_var);

    fprintf('Dormand-Prince fixed-step', slope_fix);
    fprintf('Dormand-Prince adaptive-step', slope_var);

    %% Plot
    figure();
    loglog(h_avg_fix, error_fix, 'b.', 'MarkerSize', 15, ...
        'DisplayName', 'DP Fixed Num.');
    hold on;
    loglog(h_avg_fix, yfit_fix, 'k-', 'LineWidth', 2, ...
        'DisplayName', sprintf('DP Fixed Fit', slope_fix));
    loglog(h_avg_var, error_var, 'r.', 'MarkerSize', 15, ...
        'DisplayName', 'DP Adaptive Num.');
    loglog(h_avg_var, yfit_var, 'k-', 'LineWidth', 2, ...
        'DisplayName', sprintf('DP Adaptive Fit', slope_var));
    hold off;
    grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log',"TickLabelInterpreter",'latex', 'Fontsize',12);
    title('Global Truncation Error vs. Average Step Size (Dormand-Prince)', 'FontSize', 20, 'Interpreter', 'Latex');
    xlabel('$h_{avg}(-)$', 'Interpreter', 'Latex', 'FontSize', 15);
    ylabel('$\epsilon (-)$', 'Interpreter', 'Latex', 'FontSize', 15);
    legend('Location', 'northwest', 'Interpreter', 'Latex', 'FontSize',15);

    %% Number of evaluations: fits (same error-based fit region as above)
    mask_fix_e = isfinite(error_fix) & isfinite(evals_fix) & ...
                 error_fix >= err_bounds_fixed(1) & error_fix <= err_bounds_fixed(2);
    mask_var_e = isfinite(error_var) & isfinite(evals_var) & ...
                 error_var >= err_bounds_adaptive(1) & error_var <= err_bounds_adaptive(2);

    [slope_fix_e, yfit_fix_e] = fit_loglog(evals_fix, error_fix, mask_fix_e);
    [slope_var_e, yfit_var_e] = fit_loglog(evals_var, error_var, mask_var_e);

    fprintf('Dormand-Prince fixed-step evals', slope_fix_e);
    fprintf('Dormand-Prince adaptive-step evals', slope_var_e);

    %% Plot: error vs. number of evaluations
    figure();
    loglog(evals_fix, error_fix, 'b.', 'MarkerSize', 15, ...
        'DisplayName', 'DP Fixed Num.');
    hold on;
    loglog(evals_fix, yfit_fix_e, 'k-', 'LineWidth', 2, ...
        'DisplayName', sprintf('DP Fixed Fit', slope_fix_e));
    loglog(evals_var, error_var, 'r.', 'MarkerSize', 15, ...
        'DisplayName', 'DP Adaptive Num.');
    loglog(evals_var, yfit_var_e, 'k-', 'LineWidth', 2, ...
        'DisplayName', sprintf('DP Adaptive Fit', slope_var_e));
    hold off;
    grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log',"TickLabelInterpreter",'latex', 'Fontsize',12);
    title('Global Truncation Error vs. Number of Evaluations (Dormand-Prince)', 'Interpreter', 'Latex', 'FontSize', 19);
    xlabel('$NumEvals(-)$', 'Interpreter', 'Latex', 'FontSize', 15);
    ylabel('$\epsilon (-)$', 'Interpreter', 'Latex', 'FontSize', 15);
    legend('Location', 'southwest', 'Interpreter', 'Latex', 'FontSize',15);

    %% Failure rate vs. average step size (adaptive only: fixed step never rejects)
    % A log axis can't show zero, and if EVERY failure rate is zero there is no
    % positive data at all, so MATLAB can't build log axes and falls back to
    % linear ones. To keep the axes log-log, zero-failure points are drawn as
    % hollow circles at a floor value just below the smallest nonzero rate.
    zero_mask = isfinite(fail_var) & isfinite(h_avg_var) & fail_var == 0;
    pos_mask  = isfinite(fail_var) & isfinite(h_avg_var) & fail_var > 0;
    if any(pos_mask)
        fail_floor = 0.5 * min(fail_var(pos_mask));
    else
        fail_floor = 1e-2;   % arbitrary display value when nothing ever failed
    end

    figure();
    loglog(h_avg_var(pos_mask), fail_var(pos_mask), 'r.', 'MarkerSize', 15, ...
        'DisplayName', 'DP Adaptive Num.');
    grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log',"TickLabelInterpreter",'latex', 'Fontsize',12);
    title('Failure Rate vs. Average Step Size (Dormand-Prince)', 'Interpreter', 'Latex', 'FontSize', 20);
    xlabel('$h_{avg}(-)$', 'Interpreter', 'Latex', 'FontSize', 15);
    ylabel('Failure Rate (-)', 'Interpreter', 'Latex', 'FontSize', 15);
    legend('Location', 'southeast', 'Interpreter', 'Latex', 'FontSize', 15);
end

function [slope, yfit] = fit_loglog(h, err, mask)
    slope = NaN;
    yfit = NaN(size(h));
    if sum(mask) >= 2
        c = polyfit(log10(h(mask)), log10(err(mask)), 1);
        slope = c(1);
        yfit(mask) = 10.^polyval(c, log10(h(mask)));
    else
        warning('Fewer than 2 points in the fit region; adjust the error bounds.');
    end
end
