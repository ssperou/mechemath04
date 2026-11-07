function plot_orbitdistance()
    %Initialize relevant structs
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

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

    %Set initial values
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,35,100);
    h_ref = 1e-4;
    t_span = [0 35];
    grav_func = @(t, V) gravity_rate_func(t, V, orbit_params);

    p = 5;
    h_ref = 1e-2;
    e_des = 1e-6;
    grow_lim = 2;
    sun = [0 0];

    %Analytical solution
    V_list = compute_planetary_motion(t_range, V0, orbit_params);
    X0 = V_list(:, 1);

    %Run integration
    [tlist_DP, Xlist_DP] = explicit_RK_variable_step_integration(grav_func, t_span, V0, h_ref,  DormandPrince, p, e_des, grow_lim);

    %Pull lists
    position = Xlist_DP(:, 1:2);
    h_list = diff(tlist_DP(:));
    dist = vecnorm(position - sun, 2, 2);
    dist = dist(2:end);

    figure();
    scatter(h_list, dist, 15, 'filled');
    grid on
    title('Step Size vs. Distance from Sun (Dormand-Prince)', "Interpreter", "Latex", "FontSize", 17)
    xlabel ("Step Size (-)", "Interpreter", "Latex", "FontSize", 12)
    ylabel ("Distance (-)", "Interpreter", "Latex", "FontSize", 12)
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

end