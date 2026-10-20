function error_val = kepler_func_hyperbola(y,a,b,target_area)
    c = sqrt(a^2+b^2);
    q = y/b;
    hyperbolic_angle = log(q+sqrt(1+q^2));

    error_val = y*c/2-a*b*hyperbolic_angle/2 - target_area;
end