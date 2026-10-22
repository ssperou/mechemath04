function plottingRKGlobalTruncation()
    %Initialize structs
    ForwardEuler = struct();
    ForwardEuler.A = [0, 0;2/3, 0;];
    ForwardEuler.B = [0.25, 0.75];
    ForwardEuler.C = [0, 2/3];

    ExplicitMidpoint = struct();
    ExplicitMidpoint.A = [0, 0;0.5, 0];
    ExplicitMidpoint.B = [0, 1];
    ExplicitMidpoint.C = [0, 0.5];

    HouwenWray = struct();
    HouwenWray.A = [0, 0, 0;8/15, 0, 0;1/4, 5/12, 0];
    HouwenWray.B = [0.25, 0, 0.75];
    HouwenWray.C = [0, 8/15, 2/3];

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

    X0 = solution(t_span(1), V0, orbit_params);
    X0 = X0(:);
    h_list = logspace(-3, 0, 40);
    log_h_list = log10(h_list)';

    error_FE = zeros(40, 1);
    num_iter_FE = zeros(40, 1);
    error_EM = zeros(40, 1);
    num_iter_EM = zeros(40, 1);
    error_HW = zeros(40, 1);
    num_iter_HW = zeros(40, 1);
    %Forward Euler
    for i = 1:40
        h = h_list(i);
        [~,X_list_FE,~, num_evals_FE] = explicit_RK_fixed_step_integration(gravity_rate_func_in,t_span,X0,h,ForwardEuler);
        x_real = solution(t_span(end), V0, orbit_params);
        error_FE(i) = norm(X_list_FE(end,:).' - x_real);
        num_iter_FE(i) = num_evals_FE;
    end
    log_error_FE = log10(error_FE);
    %Explicit Midpoint
    for i = 1:40
        h = h_list(i);
        [~,X_list_EM,~, num_evals_EM] = explicit_RK_fixed_step_integration(gravity_rate_func_in,t_span,X0,h,ExplicitMidpoint);
        x_real = solution(t_span(end), V0, orbit_params); x_real = x_real(:);
        error_EM(i) = norm(X_list_EM(end,:).' - x_real);
        num_iter_EM(i) = num_evals_EM;
    end
    log_error_EM = log10(error_EM);

    %Houwen/Wray
    for i = 1:40
        h = h_list(i);
        [~,X_list_HW,~, num_evals_HW] = explicit_RK_fixed_step_integration(gravity_rate_func_in,t_span,X0,h,HouwenWray);
        x_real = solution(t_span(end), V0, orbit_params);
        error_HW(i) = norm(X_list_HW(end,:).' - x_real);
        num_iter_HW(i) = num_evals_HW;
    end
    log_error_HW = log10(error_HW);

    %Step Size Global Trunction
    x_FE = polyfit(log_h_list, log_error_FE, 1);
    y_line_FE = 10.^polyval(x_FE, log_h_list);
    x_EM = polyfit(log_h_list, log_error_EM, 1);
    y_line_EM = 10.^polyval(x_EM, log_h_list);
    x_HW = polyfit(log_h_list, log_error_HW, 1);
    y_line_HW = 10.^polyval(x_HW, log_h_list);

    p_FE = x_FE(1) %p value
    p_EM = x_EM(1) %p value
    p_HW = x_HW(1) %p value
    
    %Plotting
    figure();
    loglog(h_list, y_line_FE, 'k-');
    hold on;
    loglog(h_list, error_FE, 'r.');
    loglog(h_list, y_line_EM, 'k-');
    loglog(h_list, error_EM, 'b.');
    loglog(h_list, y_line_HW, 'k-');
    loglog(h_list, error_HW, 'g.');
    hold off;
    title("Global Truncation vs. Step Size", "Interpreter", "Latex");
    xlabel('$h(-)$', 'Interpreter', 'Latex');
    ylabel('$\epsilon (-)$', 'Interpreter', 'Latex');
    legend('Ralston Fit', 'Forward Euler Num.', 'Explicit Midpoint Fit', 'Explicit Midpoint Num.', 'Houwen Wray Fit', 'Houwen Wray Num.', 'Interpreter', 'Latex', 'Location', 'southeast');

    %Num Evals
    log_evals_FE = log10(num_iter_FE);
    x_evals_FE = polyfit(log_evals_FE, log_error_FE, 1);
    y_fit_FE = 10.^polyval(x_evals_FE, log_evals_FE);
    log_evals_EM = log10(num_iter_EM);
    x_evals_EM = polyfit(log_evals_EM, log_error_EM, 1);
    y_fit_EM = 10.^polyval(x_evals_EM, log_evals_EM);
    log_evals_HW = log10(num_iter_HW);
    x_evals_HW = polyfit(log_evals_HW, log_error_HW, 1);
    y_fit_HW = 10.^polyval(x_evals_HW, log_evals_HW);

    p_FE = x_evals_FE(1) %p value
    p_EM = x_evals_EM(1) %p value
    p_HW = x_evals_HW(1) %p value

    figure();
    loglog(num_iter_FE, y_fit_FE, 'k-');
    hold on;
    loglog(num_iter_FE, error_FE, 'r.');
    loglog(num_iter_EM, y_fit_EM, 'k-');
    loglog(num_iter_EM, error_EM, 'b.');
    loglog(num_iter_HW, y_fit_HW, 'k-');
    loglog(num_iter_HW, error_HW, 'g.');
    hold off;
    title("Global Truncation vs. Num Evals", "Interpreter", "Latex");
    xlabel('$NumEvals(-)$', 'Interpreter', 'Latex');
    ylabel('$\epsilon (-)$', 'Interpreter', 'Latex');
    legend('Ralston Fit', 'Ralston Num.', 'Explicit Midpoint Fit', 'Explicit Midpoint Num.', 'Houwen Wray Fit', 'Houwen Wray Num.', 'Interpreter', 'Latex', 'Location', 'southwest');
end