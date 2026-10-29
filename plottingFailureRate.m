function plottingFailureRate()
    testfunc = ;
    sol = ;
    
    Ralston = struct();
    Ralston.A = [0,0;2/3,0];
    Ralston.B = [1/4,3/4];
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
    num_iter = 40;
    t_span = [0 2*pi];
    solution = @compute_planetary_motion;
    gravity_rate_func_in = @(t, V) gravity_rate_func(t, V, orbit_params);

    X0 = V0;
    X0 = X0(:)';
    h_list = logspace(-3, 0, num_iter);
    log_h_list = log10(h_list)';
    e_des = 9;

    evals_list = zeros();
    failed_list = zeros();

    for i = 1:num_iter
        [t_list,X_list,h_avg, num_evals, failed] = explicit_RK_variable_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct,p,error_desired)


end