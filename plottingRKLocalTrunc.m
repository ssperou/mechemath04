% plot local truncation errors for explicit runge kutta methods
function plottingRKLocalTrunc(sol, t_ref, num_iter)
    arguments
        sol (1,:) function_handle = @compute_planetary_motion;
        t_ref (1,:) double = 0.439
        num_iter (1,:) double = 100
    end
    %Input necessary structs
    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
                       1/5, 0, 0, 0,0,0,0;...
                       3/40, 9/40, 0, 0, 0, 0,0;...
                       44/45, -56/15, 32/9, 0, 0, 0,0;...
                       19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
                       9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
                       35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    %Initialize values
    func = @(t, V) gravity_rate_func(t, V, orbit_params);
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0;y0;dxdt0;dydt0];

    XA = sol(t_ref,V0,orbit_params);
    XA = XA(:);
    h_list = logspace(-3,1,num_iter);

    embedded_error = zeros(num_iter,1);
    xb1_error = zeros(num_iter,1);
    xb2_error = zeros(num_iter,1);
   

    for i = 1:num_iter
        h = h_list(i);
        [XB1, XB2, num_evals] = RK_step_embedded(func, t_ref, XA, h, DormandPrince);

        x_true = sol(t_ref + h, V0, orbit_params);

        embedded_error(i) = norm(XB1 - XB2);

        xb1_error(i) = norm(XB1 - x_true);
        xb2_error(i) = norm(XB2 - x_true);
       
    end

    % %%%%%%% P-VALUE CALCULATION
    % h_filter = h_list <= 1e-2;
    % log_h_list = log10(h_list(h_filter))';
    % 
    % % solution p-value calculation
    % log_solution_list = log10(solution_list(h_filter));
    % solution_coeff = polyfit(log_h_list, log_solution_list, 1);
    % 
    % % Explicit Midpoint p-value calculation
    % log_explicit_midpoint_error = log10(explicit_midpoint_error(h_filter));
    % explicit_midpoint_coeff = polyfit(log_h_list, log_explicit_midpoint_error, 1);
    % 
    % % Forward Euler p-value calculation
    % log_forward_euler_error = log10(forward_euler_error(h_filter));
    % forward_euler_coeff = polyfit(log_h_list, log_forward_euler_error, 1);
    % 
    % % % % Implicit Midpoint p-value calculation
    % log_implicit_midpoint_error = log10(implicit_midpoint_error(h_filter));
    % implicit_midpoint_coeff = polyfit(log_h_list, log_implicit_midpoint_error, 1);
    % 
    % % Backward Euler p-value calculation
    % log_backward_euler_error = log10(backward_euler_error(h_filter));
    % backward_euler_coeff = polyfit(log_h_list, log_backward_euler_error, 1);
    % 
    % % % P-VALUES
    % solution_p = solution_coeff(1)
    % forward_p = forward_euler_coeff(1)
    % backward_p = backward_euler_coeff(1)
    % explicit_p = explicit_midpoint_coeff(1)
    % implicit_p = implicit_midpoint_coeff(1)
    % 
    % 
    % 
    % % FIT LINE CODE
    % explicit_midpoint_fit = 10.^polyval(explicit_midpoint_coeff, log10(h_list));
    % forward_euler_fit = 10.^polyval(forward_euler_coeff, log10(h_list));
    % solution_fit = 10.^polyval(solution_coeff, log10(h_list));
    % 
    % %%%%%% FIT LINE PLOTS
    % loglog(10.^log10(h_list), explicit_midpoint_fit, '-', 'Linewidth', 5, 'Color', [0.5 1 0.5], 'Displayname',"Exp. Mid. Fit");
    % hold on;
    % grid on;
    % loglog(10.^log10(h_list), forward_euler_fit, '-', 'Linewidth', 5, 'Color', [1 0.5 0.5], 'Displayname',"Forw. Euler Fit");
    % loglog(10.^log10(h_list), solution_fit, '-', 'Linewidth', 5, 'Color', [0.6 0.6 1], 'Displayname',"Abs. Diff. Fit");


    %%%%%%% STEP SIZE PLOTS
    figure;
    loglog(h_list, embedded_error, 'k-', 'markersize', 25, 'Displayname',"|XB1-XB2|");
    hold on;
    loglog(h_list, xb1_error, 'b.', 'markersize', 25, 'Displayname',"XB1");
    loglog(h_list, xb2_error, 'g.', 'markersize', 25, 'Displayname',"XB2");
   
    grid on;

    title('Local Truncation Error $\epsilon$ for Four Methods vs Timesteps', 'Interpreter', 'Latex', 'FontSize',16);
    legend('Location','southeast', 'Interpreter', 'Latex','FontSize',12);
    xlabel("Step size $h$ (-)", 'Interpreter', 'Latex','FontSize',16);
    ylabel("Local Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
    set(gca,"TickLabelInterpreter",'latex');
    %fontsize(scale=2)  
   
end