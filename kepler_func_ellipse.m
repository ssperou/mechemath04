function error_val = kepler_func_ellipse(theta,a,b,target_area)
    i = sqrt(-1);

    c = sqrt(a^2-b^2);
    r = b^2/(a-c*cos(theta));

    x = (r*cos(theta)-c)/a;
    y = (r*sin(theta))/b;
    psi = angle(x+i*y);

    error_val = a*b*psi/2 + (c/2)*r*sin(theta) - target_area;
end
