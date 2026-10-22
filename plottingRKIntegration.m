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

    %  Ralston's method
    Ralston = struct();
    Ralston.A = [0,0;2/3,0];
    Ralston.B = [1/4,3/4];
    Ralston.C = [0, 2/3];

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

    t_range = linspace(0,36,100);
    t_span = [0 36];
    h_ref = 0.03;
    V_list = compute_planetary_motion(t_range, V0, orbit_params);
    X0 = V_list(1, :);
    gravity_rate_func_in = @(t, V) gravity_rate_func(t, V, orbit_params);
    
    %Ralston Euler Integration
    [t_list_FE,X_list_FE,h_avg_FE, num_evals_FE] = explicit_RK_fixed_step_integration(gravity_rate_func_in,t_span,X0,h_ref,Ralston);
    
    %Explicit Midpoint
    [t_list_EM,X_list_EM,h_avg_EM, num_evals_EM] = explicit_RK_fixed_step_integration(gravity_rate_func_in,t_span,X0,h_ref,ExplicitMidpoint)
    
    %Van der Houwen/Wray
    [t_list_HW,X_list_HW,h_avg_HW, num_evals_HW] = explicit_RK_fixed_step_integration(gravity_rate_func_in,t_span,X0,h_ref,HouwenWray)
    



    fig1 = figure(1);
    tiledlayout(2,1,"TileSpacing",'compact','Padding','compact');

    nexttile();

    hold on;

    %Plotting integration
    h1 = plot(t_range, V_list(:, 1), 'k-', 'LineWidth', 4); %Analytical Solution
    hold on;
    h2 = plot(t_list_FE, X_list_FE(:, 1), 'r-', 'LineWidth', 2.75);
    h3 = plot(t_list_EM, X_list_EM(:,1), 'b--', 'LineWidth', 2.5);
    h4 = plot(t_list_HW, X_list_HW(:, 1), 'g--', 'LineWidth', 2.25);
    title("Runge Kutta Integration (X Component)", 'Interpreter', 'Latex','FontSize',18);
    xlabel("$Time\ (-)$", 'Interpreter', 'Latex','FontSize',16); xlabel("$Time\ (-)$", 'Interpreter', 'Latex','FontSize',16);
    ylabel("$X(t)$", 'Interpreter', 'Latex','FontSize',16);
    set(gca,"TickLabelInterpreter",'latex', 'FontSize',13);
    axis([0, 36, -5, 10])
    %legend([h1, h2, h3, h4], "Analytical", "Ralston", "Explicit Midpoint", "Houwen/Wray", "Interpreter", "Latex", "Location", "northeast");
    nexttile();
    hold on;

    %Plotting y integration
    h1 = plot(t_range, V_list(:, 2), 'k-','LineWidth', 4); %Analytical Solution
    hold on;
    h2 = plot(t_list_FE, X_list_FE(:, 2), 'r-', 'LineWidth', 2.75);
    h3 = plot(t_list_EM, X_list_EM(:,2), 'b--', 'LineWidth', 2.5);
    h4 = plot(t_list_HW, X_list_HW(:, 2), 'g--', 'LineWidth', 2.25);
    title("Runge Kutta Integration (Y Component)", 'Interpreter', 'Latex','FontSize',18);
    xlabel("$Time\ (-)$", 'Interpreter', 'Latex','FontSize',16); xlabel("$Time\ (-)$", 'Interpreter', 'Latex','FontSize',16);
    set(gca,"TickLabelInterpreter",'latex','FontSize',13);
    ylabel("$X(t)$", 'Interpreter', 'Latex','FontSize',16);
    axis([0, 36, -5, 5])
 
    
    my_legend = legend([h1, h2, h3, h4], "Analytical", "Ralston", "Explicit Midpoint", "Houwen/Wray", "Interpreter", "Latex", "Location", "northeast");
    my_legend.Layout.Tile = 'south';
  
    set(my_legend, 'orientation','horizontal', 'FontSize', 15)

    % figure();
    % axis equal;
    % h1 = plot(V_list(:,1), V_list(:, 2), 'k-'); %Analytical Solution
    % hold on;
    % plot(0, 0, 'ko', 'MarkerSize', 15);
    % h2 = plot(X_list_FE(:, 1), X_list_FE(:, 2), 'r-');
    % h3 = plot(X_list_EM(:, 1), X_list_EM(:,2), 'b--');
    % h4 = plot(X_list_HW(:, 1), X_list_HW(:, 2), 'g:');
    % title("Runge Kutta Integration (Y Values)", 'Interpreter', 'Latex');
    % xlabel("$t(-)$", 'Interpreter', 'Latex'); xlabel("$t(-)$", 'Interpreter', 'Latex');
    % legend([h1, h2, h3, h4], "Analytical", "Ralston Euler", "Explicit Midpoint", "Houwen/Wray", "Interpreter", "Latex", "Location", "northeast");
    % 
    % hold off;

    
end