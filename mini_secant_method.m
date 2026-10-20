function x_out = mini_secant_method(func_in,xa)
    my_tol = 1e-13;

    xb = xa+1e-2;

    ya = func_in(xa);
    yb = func_in(xb);

    count = 0;
    while count<20 && abs(yb)>my_tol
        count = count+1;

        xc = (xa*yb-xb*ya)/(yb-ya);
        yc = func_in(xc);

        xa = xb;
        ya = yb;

        xb = xc;
        yb = yc;
    end

    x_out = xb;
end