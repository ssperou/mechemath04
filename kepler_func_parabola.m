function error_val = kepler_func_parabola(y,d,target_area)
    error_val = y^3/(12*d)+y*d/4 - target_area;
end
