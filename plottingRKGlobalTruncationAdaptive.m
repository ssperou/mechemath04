function plottingRKGlobalTruncationAdaptive()
    %Initialize structs
    % ForwardEuler = struct();
    % ForwardEuler.A = [0, 0;1, 0];
    % ForwardEuler.B = [1, 0; 0.5, 0.5];
    % ForwardEuler.C = [0, 1];
    % 
    % ExplicitMidpoint = struct();
    % ExplicitMidpoint.A = [0, 0;0.5, 0];
    % ExplicitMidpoint.B = [0, 1; 1, 0];
    % ExplicitMidpoint.C = [0, 0.5];

    DormandPrince.A = [0,0,0,0,0,0,0;
        1/5, 0, 0, 0,0,0,0;...
        3/40, 9/40, 0, 0, 0, 0,0;...
        44/45, -56/15, 32/9, 0, 0, 0,0;...
        19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
        9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
        35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];
    HouwenWray = struct();
    HouwenWray.A = [0, 0, 0;8/15, 0, 0;1/4, 5/12, 0];
    HouwenWray.B = [0.25, 0, 0.75; 1/16, 15/16, 0];
    HouwenWray.C = [0, 8/15, 2/3];

    Ralston = struct();
    Ralston.A = [0,0;2/3,0];
    Ralston.B = [1/4,3/4; 1,0];
    Ralston.C = [0, 2/3];

    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    %Initiate values
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];
    t_span = [0 2*pi];
    solution = @compute_planetary_motion;
    gravity_rate_func_in = @(t, V) gravity_rate_func(t, V, orbit_params);

    X0 = V0(:);
    x_real = solution(t_span(end), V0, orbit_params);
    x_real = x_real(:);

    %Sweep settings: tolerance is swept, h is only the starting guess
    n_pts = 50;
    h_ref = 1e-2;
    growth_limit = 8;
    e_list = logspace(-3, 0, n_pts);
    e_list = e_list(end:-1:1);

    %Bounds of the linear region used for the fits (tune these by looking at the plots)
    %[min max] of h_avg for the step size and failure fits
    h_bounds_FE = [1e-4, 1e-2];
    h_bounds_EM = [1e-3, 1e-1];
    h_bounds_HW = [1e-2, 1e0];
    h_bounds_R  = [1e-3, 1e-1];
    %[min max] of function evals for the efficiency fits
    evals_bounds_FE = [100, 500];
    evals_bounds_EM = [1e2, 1e4];
    evals_bounds_HW = [1e2, 1e3];
    evals_bounds_R  = [1e2, 1e4];

    %Create lists
    error_FE = zeros(n_pts, 1);
    num_evals_FE = zeros(n_pts, 1);
    h_avg_FE = zeros(n_pts, 1);
    fail_FE = zeros(n_pts, 1);
    error_EM = zeros(n_pts, 1);
    num_evals_EM = zeros(n_pts, 1);
    h_avg_EM = zeros(n_pts, 1);
    fail_EM = zeros(n_pts, 1);
    error_HW = zeros(n_pts, 1);
    num_evals_HW = zeros(n_pts, 1);
    h_avg_HW = zeros(n_pts, 1);
    fail_HW = zeros(n_pts, 1);
    error_R = zeros(n_pts, 1);
    num_evals_R = zeros(n_pts, 1);
    h_avg_R = zeros(n_pts, 1);
    fail_R = zeros(n_pts, 1);

    %Forward Euler
    % p = 1;
    % e_list_FE = logspace(-2, 0, n_pts);
    % e_list_FE = e_list_FE(end:-1:1);
    % for i = 1:n_pts
    %     e_des = e_list_FE(i);
    %     [~,X_list_FE,h_avg, evals, failed] = explicit_RK_variable_step_integration(gravity_rate_func_in,t_span,X0,h_ref,ForwardEuler, p, e_des, growth_limit);
    %     error_FE(i) = norm(X_list_FE(end,:).' - x_real);
    %     num_evals_FE(i) = evals;
    %     h_avg_FE(i) = h_avg;
    %     fail_FE(i) = failed;
    % end

    %Explicit Midpoint
    p = 2;
    for i = 1:n_pts
        e_des = e_list(i);
        [~,X_list_EM,h_avg, evals, failed] = explicit_RK_variable_step_integration(gravity_rate_func_in,t_span,X0,h_ref,ExplicitMidpoint, p, e_des, growth_limit);
        error_EM(i) = norm(X_list_EM(end,:).' - x_real);
        num_evals_EM(i) = evals;
        h_avg_EM(i) = h_avg;
        fail_EM(i) = failed;
    end

    %Houwen/Wray
    p = 3;
    for i = 1:n_pts
        e_des = e_list(i);
        [~,X_list_HW,h_avg, evals, failed] = explicit_RK_variable_step_integration(gravity_rate_func_in,t_span,X0,h_ref,HouwenWray, p, e_des, growth_limit);
        error_HW(i) = norm(X_list_HW(end,:).' - x_real);
        num_evals_HW(i) = evals;
        h_avg_HW(i) = h_avg;
        fail_HW(i) = failed;
    end

    %Ralston's
    p = 2;
    for i = 1:n_pts
        e_des = e_list(i);
        [~,X_list_R,h_avg, evals, failed] = explicit_RK_variable_step_integration(gravity_rate_func_in,t_span,X0,h_ref,Ralston, p, e_des, growth_limit);
        error_R(i) = norm(X_list_R(end,:).' - x_real);
        num_evals_R(i) = evals;
        h_avg_R(i) = h_avg;
        fail_R(i) = failed;
    end

    %% Step Size Global Truncation (fit only the linear region)
    %filt_FE = h_avg_FE >= h_bounds_FE(1) & h_avg_FE <= h_bounds_FE(2);
    %x_FE = polyfit(log10(h_avg_FE(filt_FE)), log10(error_FE(filt_FE)), 1);
    %y_line_FE = 10.^polyval(x_FE, log10(h_avg_FE));
    filt_EM = h_avg_EM >= h_bounds_EM(1) & h_avg_EM <= h_bounds_EM(2);
    x_EM = polyfit(log10(h_avg_EM(filt_EM)), log10(error_EM(filt_EM)), 1);
    y_line_EM = 10.^polyval(x_EM, log10(h_avg_EM));
    filt_HW = h_avg_HW >= h_bounds_HW(1) & h_avg_HW <= h_bounds_HW(2);
    x_HW = polyfit(log10(h_avg_HW(filt_HW)), log10(error_HW(filt_HW)), 1);
    y_line_HW = 10.^polyval(x_HW, log10(h_avg_HW));
    filt_R = h_avg_R >= h_bounds_R(1) & h_avg_R <= h_bounds_R(2);
    x_R = polyfit(log10(h_avg_R(filt_R)), log10(error_R(filt_R)), 1);
    y_line_R = 10.^polyval(x_R, log10(h_avg_R));

    %p_FE = x_FE(1) %p value
    p_EM = x_EM(1) %p value
    p_HW = x_HW(1) %p value
    p_R = x_R(1)

    %Plotting
    figure();
    %loglog(h_avg_FE, y_line_FE, 'k-', 'DisplayName', 'Forward Euler Fit');
    %hold on;
    %loglog(h_avg_FE, error_FE, 'r.', 'MarkerSize', 15, 'DisplayName', 'Forward Euler Num.');
    loglog(h_avg_EM, y_line_EM, 'k-', 'DisplayName', 'Explicit Midpoint Fit');
    hold on;
    loglog(h_avg_EM, error_EM, 'b.', 'MarkerSize', 15, 'DisplayName', 'Explicit Midpoint Num.');
    loglog(h_avg_HW, y_line_HW, 'k-', 'DisplayName', 'Houwen Wray Fit');
    loglog(h_avg_HW, error_HW, 'g.', 'MarkerSize', 15, 'DisplayName', 'Houwen Wray Num.');
    loglog(h_avg_R, y_line_R, 'k-', 'DisplayName', 'Ralston Fit');
    loglog(h_avg_R, error_R, 'm.', 'MarkerSize', 15, 'DisplayName', 'Ralston Num.');
    hold off;
    grid on;
    title("Global Truncation vs. Average Step Size", "Interpreter", "Latex");
    xlabel('$h_{avg}(-)$', 'Interpreter', 'Latex');
    ylabel('$\epsilon (-)$', 'Interpreter', 'Latex');
    legend('Interpreter', 'Latex', 'Location', 'northwest');

    %% Num Evals (fit only the linear region)
    %filt_FE = num_evals_FE >= evals_bounds_FE(1) & num_evals_FE <= evals_bounds_FE(2);
    %x_evals_FE = polyfit(log10(num_evals_FE(filt_FE)), log10(error_FE(filt_FE)), 1);
   % y_fit_FE = 10.^polyval(x_evals_FE, log10(num_evals_FE));
    filt_EM = num_evals_EM >= evals_bounds_EM(1) & num_evals_EM <= evals_bounds_EM(2);
    x_evals_EM = polyfit(log10(num_evals_EM(filt_EM)), log10(error_EM(filt_EM)), 1);
    y_fit_EM = 10.^polyval(x_evals_EM, log10(num_evals_EM));
    filt_HW = num_evals_HW >= evals_bounds_HW(1) & num_evals_HW <= evals_bounds_HW(2);
    x_evals_HW = polyfit(log10(num_evals_HW(filt_HW)), log10(error_HW(filt_HW)), 1);
    y_fit_HW = 10.^polyval(x_evals_HW, log10(num_evals_HW));
    filt_R = num_evals_R >= evals_bounds_R(1) & num_evals_R <= evals_bounds_R(2);
    x_evals_R = polyfit(log10(num_evals_R(filt_R)), log10(error_R(filt_R)), 1);
    y_fit_R = 10.^polyval(x_evals_R, log10(num_evals_R));

    %p_FE = x_evals_FE(1) %p value
    p_EM = x_evals_EM(1) %p value
    p_HW = x_evals_HW(1) %p value
    p_R = x_evals_R(1)

    figure();
    %loglog(num_evals_FE, y_fit_FE, 'k-','DisplayName', 'Forward Euler Fit');
    %hold on;
    %loglog(num_evals_FE, error_FE, 'r.', 'MarkerSize', 15, 'DisplayName', 'Forward Euler Num.');
    loglog(num_evals_EM, y_fit_EM, 'k-', 'DisplayName', 'Explicit Midpoint Fit');
    hold on;
    loglog(num_evals_EM, error_EM, 'b.', 'MarkerSize', 15, 'DisplayName', 'Explicit Midpoint Num.');
    loglog(num_evals_HW, y_fit_HW, 'k-', 'DisplayName', 'Houwen Wray Fit');
    loglog(num_evals_HW, error_HW, 'g.', 'MarkerSize', 15, 'DisplayName', 'Houwen Wray Num.');
    loglog(num_evals_R, y_fit_R, 'k-', 'DisplayName', 'Ralston Fit');
    loglog(num_evals_R, error_R, 'm.', 'MarkerSize', 15, 'DisplayName', 'Ralston Num.');
    hold off;
    grid on;
    title("Global Truncation vs. Num Evals", "Interpreter", "Latex");
    xlabel('$NumEvals(-)$', 'Interpreter', 'Latex');
    ylabel('$\epsilon (-)$', 'Interpreter', 'Latex');
    legend('Interpreter', 'Latex', 'Location', 'southwest');

    %% Failure Rate vs Step Size (zero failure rates can't be logged, so drop them)
  %% Failure Rate vs Step Size
    % filt_FE = h_avg_FE >= h_bounds_FE(1) & h_avg_FE <= h_bounds_FE(2) & fail_FE > 0;
    % y_fail_FE = NaN(size(h_avg_FE));
    % if sum(filt_FE) >= 3
    %     x_fail_FE = polyfit(log10(h_avg_FE(filt_FE)), log10(fail_FE(filt_FE)), 1);
    %     y_fail_FE = 10.^polyval(x_fail_FE, log10(h_avg_FE));
    % end

    filt_EM = h_avg_EM >= h_bounds_EM(1) & h_avg_EM <= h_bounds_EM(2) & fail_EM > 0;
    y_fail_EM = NaN(size(h_avg_EM));
    if sum(filt_EM) >= 3
        x_fail_EM = polyfit(log10(h_avg_EM(filt_EM)), log10(fail_EM(filt_EM)), 1);
        y_fail_EM = 10.^polyval(x_fail_EM, log10(h_avg_EM));
    end

    filt_R = h_avg_R >= h_bounds_R(1) & h_avg_R <= h_bounds_R(2) & fail_R > 0;
    y_fail_R = NaN(size(h_avg_R));
    if sum(filt_R) >= 3
        x_fail_R = polyfit(log10(h_avg_R(filt_R)), log10(fail_R(filt_R)), 1);
        y_fail_R = 10.^polyval(x_fail_R, log10(h_avg_R));
    end

    filt_HW = h_avg_HW >= h_bounds_HW(1) & h_avg_HW <= h_bounds_HW(2) & fail_HW > 0;
    y_fail_HW = NaN(size(h_avg_HW));
    if sum(filt_HW) >= 3
        x_fail_HW = polyfit(log10(h_avg_HW(filt_HW)), log10(fail_HW(filt_HW)), 1);
        y_fail_HW = 10.^polyval(x_fail_HW, log10(h_avg_HW));
    end

    figure();
    %loglog(h_avg_FE, y_fail_FE, 'k-', 'DisplayName', 'Forward Euler Fit');
    %hold on;
    %loglog(h_avg_FE, fail_FE, 'r.', 'MarkerSize', 15, 'DisplayName', 'Forward Euler Num.');
    loglog(h_avg_EM, y_fail_EM, 'k-', 'DisplayName', 'Explicit Midpoint Fit');
    hold on;
    loglog(h_avg_EM, fail_EM, 'b.', 'MarkerSize', 15, 'DisplayName', 'Explicit Midpoint Num.');
    loglog(h_avg_HW, y_fail_HW, 'k-', 'DisplayName', 'Houwen Wray Fit');
    loglog(h_avg_HW, fail_HW, 'g.', 'MarkerSize', 15, 'DisplayName', 'Houwen Wray Num.');
    loglog(h_avg_R, y_fail_R, 'k-', ...
        'DisplayName', 'Ralston Fit');
    loglog(h_avg_R, fail_R, 'm.', 'MarkerSize', 15, 'DisplayName', 'Ralston Num.');
    hold off;
    grid on;
    title("Failure Rate vs. Average Step Size", "Interpreter", "Latex");
    xlabel('$h_{avg}(-)$', 'Interpreter', 'Latex');
    ylabel('Failure Rate (-)', 'Interpreter', 'Latex');
    legend('Interpreter', 'Latex', 'Location', 'northwest');
end