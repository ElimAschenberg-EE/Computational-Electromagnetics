function Q = Q_function(I, theta)

    global N d

    Q = 0;

    for n = 2:N

        Q = Q + I(n - 1) * psi_funct(3, n, theta, 0) + (I(n) - I(n - 1)) / d * psi_funct(4, n, theta, 0);

    end

end