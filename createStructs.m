function createStructs();
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
end