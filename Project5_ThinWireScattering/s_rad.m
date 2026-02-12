function s_rad = s_rad_function(I, theta)

    global w mu0 E0

    s_rad = (w * mu0 * sin(theta) * abs(Q(I, theta)))^2 / (4 * pi * E0^2);

end