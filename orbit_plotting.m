function orbit_plotting()
    %establish relevant structs
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

    HeunEuler = struct();
    HeunEuler.C = [0,1];
    HeunEuler.B = [1/2,1/2;1,0];
    HeunEuler.A = [0,0;1,0];

    Bogacki = struct();
    Bogacki.C = [0,1/2, 3/4, 1];
    Bogacki.B = [2/9, 1/3, 4/9, 0; 7/24, 1/4, 1/3, 1/8];
    Bogacki.A = [0,0,0,0; 1/2,0,0,0; 0,3/4,0,0; 2/9,1/3, 4/9, 0];

    %Find initial values
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,35,100);
    h_ref = 1e-4;
    t_span = [0 35];
    grav_func = @(t, V) gravity_rate_func(t, V, orbit_params);

    %Analytical solution
    V_list = compute_planetary_motion(t_range, V0, orbit_params);
    X0 = V_list(:, 1);

    %Run each method's integration
    [tlist_DP, Xlist_DP] = explicit_RK_fixed_step_integration(grav_func, t_span, V0, h_ref, DormandPrince);
    [tlist_HE, Xlist_HE] = explicit_RK_fixed_step_integration(grav_func, t_span, V0, h_ref, HeunEuler);
    [tlist_Bo, Xlist_Bo] = explicit_RK_fixed_step_integration(grav_func, t_span, V0, h_ref, Bogacki);

    %Plot time vs position
    figure(1)
    plot(t_range, V_list, 'k--');
    hold on
    plot(tlist_DP, Xlist_DP, 'r--');
    plot(tlist_HE, Xlist_HE, 'b--');
    plot(tlist_Bo, Xlist_Bo, 'g--');
    hold off

    %Plot orbit path
    figure(2)
    axis equal; axis square;
    axis([-20,20,-20,20]);
    hold ;
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(V_list(:, 1), V_list(:, 2), 'k--');
    plot(Xlist_DP(:, 1), Xlist_DP(:, 2), 'r--');
    plot(Xlist_HE(:, 1), Xlist_HE(:, 2), 'b--');
    plot(Xlist_Bo(:, 1), Xlist_Bo(:, 2), 'g--');
    hold off

end