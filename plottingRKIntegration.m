function plottingRKIntegration()
    %Initialize structs
    ForwardEuler = struct();
    ForwardEuler.A = [0];
    ForwardEuler.B = [1];
    ForwardEuler.C = [0];

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

    %Initialize values
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];

    t_range = linspace(0,30,100);
    t_span = [0 30];
    h_ref = 0.3;
    V_list = compute_planetary_motion(t_range, V0, orbit_params);
    X0 = V_list(1, :);
    gravity_rate_func_in = @(t, V) gravity_rate_func(t, V, orbit_params);

    %Forward Euler Integration
    [t_list_FE,X_list_FE,h_avg_FE, num_evals_FE] = explicit_RK_fixed_step_integration(gravity_rate_func_in,t_span,X0,h_ref,ForwardEuler);
    
    %Explicit Midpoint
    [t_list_EM,X_list_EM,h_avg_EM, num_evals_EM] = explicit_RK_fixed_step_integration(gravity_rate_func_in,t_span,X0,h_ref,ExplicitMidpoint);
    
    %Van der Houwen/Wray
    [t_list_HW,X_list_HW,h_avg_HW, num_evals_HW] = explicit_RK_fixed_step_integration(gravity_rate_func_in,t_span,X0,h_ref,HouwenWray)

    %Plotting integration
    figure();
    h1 = plot(t_range, V_list(:, 1), 'k-'); %Analytical Solution
    hold on;
    h2 = plot(t_list_FE, X_list_FE(:, 1), 'r--');
    h3 = plot(t_list_EM, X_list_EM(:,1), 'b--');
    h4 = plot(t_list_HW, X_list_HW(:, 1), 'g--');
    title("Runge Kutta Integration (X Values)", 'Interpreter', 'Latex');
    xlabel("$t(-)$", 'Interpreter', 'Latex'); xlabel("$t(-)$", 'Interpreter', 'Latex');
    legend([h1, h2, h3, h4], "Analytical", "Forward Euler", "Explicit Midpoint", "Houwen/Wray", "Interpreter", "Latex", "Location", "northeast");
    hold off;

    %Plotting orbit?
    figure();
    h1 = plot(t_range, V_list(:, 2), 'k-'); %Analytical Solution
    hold on;
    h2 = plot(t_list_FE, X_list_FE(:, 2), 'r--');
    h3 = plot(t_list_EM, X_list_EM(:,2), 'b--');
    h4 = plot(t_list_HW, X_list_HW(:, 2), 'g--');
    title("Runge Kutta Integration (Y Values)", 'Interpreter', 'Latex');
    xlabel("$t(-)$", 'Interpreter', 'Latex'); xlabel("$t(-)$", 'Interpreter', 'Latex');
    legend([h1, h2, h3, h4], "Analytical", "Forward Euler", "Explicit Midpoint", "Houwen/Wray", "Interpreter", "Latex", "Location", "northeast");
    hold off;
end