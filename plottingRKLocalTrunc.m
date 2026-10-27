% plot local truncation errors for explicit runge kutta methods
function plottingRKLocalTrunc(sol, t_ref, num_iter)
    arguments
        sol (1,:) function_handle = @compute_planetary_motion;
        t_ref (1,:) double = 0.439
        num_iter (1,:) double = 200
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
    h_list = logspace(-4,0,num_iter);

    embedded_error = zeros(num_iter,1);
    xb1_error = zeros(num_iter,1);
    xb2_error = zeros(num_iter,1);
    diff_error = zeros(num_iter,1);
   

    for i = 1:num_iter
        h = h_list(i);
        [XB1, XB2, num_evals] = RK_step_embedded(func, t_ref, XA, h, DormandPrince);

        x_true = sol(t_ref + h, V0, orbit_params);

        embedded_error(i) = norm(XB1 - XB2);

        xb1_error(i) = norm(XB1 - x_true);
        xb2_error(i) = norm(XB2 - x_true);
        diff_error(i) = norm(XB1 - XB2);

        timestep_error(i) = norm(x_true - XA);
       
    end

    %%%%%%% P-VALUE CALCULATIONS
    h_filter = (h_list >= 1e-1) & (h_list <= 1e0);
    log_h_list = log10(h_list(h_filter))';

    % XB1 error p-value calculation
    log_xb1_error = log10(xb1_error(h_filter));
    XB1_coeff = polyfit(log_h_list, log_xb1_error,1);
    XB1_p = XB1_coeff(1)

    % XB2 error p-value calculation
    log_xb2_error = log10(xb2_error(h_filter));
    XB2_coeff = polyfit(log_h_list, log_xb2_error,1);
    XB2_p = XB2_coeff(1)

    % Embedded diff error p-value calculation
    log_diff_error = log10(diff_error(h_filter));
    diff_coeff = polyfit(log_h_list, log_diff_error,1);
    diff_p = diff_coeff(1)


    %%%%%%% STEP SIZE PLOTS
    figure(1)
    loglog(h_list, xb1_error, 'b.', 'markersize', 15, 'Displayname',"$X_{B1}$");
    hold on;
    loglog(h_list, xb2_error, 'g.', 'markersize', 15, 'Displayname',"$X_{B2}$");
    loglog(h_list, embedded_error, 'r--', 'linewidth', 2,'Displayname',"|$X_{B1}-X_{B2}$|");
    loglog(h_list, timestep_error, 'k--', 'linewidth',2, 'Displayname',"$|f(t_{\rm ref}+h)-f(t_{\rm ref})|$");
   
    grid on;
    axis([1e-4 1e0 1e-18 1e3])

    legend('Location','southeast', 'Interpreter', 'Latex','FontSize',16);
    set(gca,"TickLabelInterpreter",'latex', 'fontsize', 14);
    title('Local Truncation Error $\epsilon$ for Embedded RK Method vs. Step Size: Dormand Prince', 'Interpreter', 'Latex', 'FontSize',20);
    xlabel("Step size $h$ (-)", 'Interpreter', 'Latex','FontSize',18);
    ylabel("Local Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',18);
    %fontsize(scale=2)  
   
end