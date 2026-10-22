% integration plots for provided method ("Explicit Midpoint",
% "Forward Euler", "Backward Euler", "Implicit Midpoint")
% orbit_params = struct();
    % orbit_params.m_sun = 1;
    % orbit_params.m_planet = 1;
    % orbit_params.G = 40;
    % x0 = 8;
    % y0 = 0;
    % dxdt0 = 0;
    % dydt0 = 1.5;
    % 
    % V0 = [x0;y0;dxdt0;dydt0];
    % t_range = linspace(0,30,100);
    % V_list = compute_planetary_motion(t_range,V0,orbit_params);
function [t_list, X_list] = plot_fixed_integration(method, href)
    arguments
        method (1,:) string = 'Explicit Midpoint'
        href (1,:) double = [0.45 2 3]
    end
    tspan = [0, 10];
    X0 = 1;
    rate_func = @gravity_rate_func;
    func_sol = @compute_planetary_motion;
    colors = [0 0 0.7; 0.5 0.25 0; 0.75 0.6 0];

    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;



    ExplicitMidpoint = struct();
    ExplicitMidpoint.A = [0, 0;0.5, 0];
    ExplicitMidpoint.B = [0, 1];
    ExplicitMidpoint.C = [0, 0.5];

    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;

    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,30,100);

    clf
    figure(1);
    hold on;
    grid on;

    for i = 1:length(href)
        h_ref = href(i);
        if method == "Explicit Midpoint"
            [t_list,X_list] = explicit_RK_fixed_step_integration( ...
                rate_func, tspan, X0, h_ref, ExplicitMidpoint);
        elseif method == "Forward Euler"
            [t_list,X_list] = forward_euler_fixed_step_integration( ...
                rate_func, tspan, X0, h_ref);
        else
            disp("not a valid method. try again")
        end 
     
     
    plot(t_list, X_list, '-o','LineWidth', 2,'MarkerFaceColor', 'k', ... 
    'MarkerSize', 3, 'MarkerEdgeColor', colors(i,:),'DisplayName', sprintf('h = %g ', h_ref));
        
    end
    
    %plot 
    xlabel('Time (-)', 'Interpreter', 'Latex','FontSize',32);
    ylabel('$X(t)$ (-)', 'Interpreter', 'Latex','FontSize',32);
    hold on;
    
    sol_t = t_range;
    x_sols = func_sol(sol_t,V0,orbit_params);

    psol = plot(sol_t(:,2), x_sols(:,2), 'LineWidth', 3, 'DisplayName', 'Analytical Solution');
    uistack(psol, 'bottom');
    legend('Location', 'southwest','Interpreter', 'Latex','FontSize',16);
    hold off;
    title(sprintf('%s Approximation vs Solution', method), 'Interpreter', 'Latex','FontSize',34);
    set(gca,"TickLabelInterpreter",'latex','FontSize',18);
    fontsize(scale=1.75)

end