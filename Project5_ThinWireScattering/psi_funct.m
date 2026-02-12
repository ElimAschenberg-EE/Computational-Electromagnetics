function val = psi_funct(type, n, theta, m)

    global epsilon0 mu0 w beta d j a

    switch type

        case 1

            if m ~= n         % ~= means logic NOT
                
                R = @(x) sqrt(( (n - m) * d - x ).^2 + a.^2);
                psi_1_function = @(x) (1 - abs(x) / d) .* (exp (-j * beta * R(x)) ./ R(x));

                val = integral(psi_1_function, -d, d);

            else

                R = @(x) sqrt(x.^2 + a^2);
                psi_1_function1 = @(x) (1 - x ./ d) .* exp(-j * beta * R(x)) ./ R(x) - cos(beta * a) ./ R(x);

                val = 2 * integral(psi_1_function1, 0, d) + 2 * cos(beta * a) * log ((d + sqrt(d^2 + a^2)) / a);
            
            end

        case 2

            R0 = sqrt(((n - m) * d)^2 + a^2);
            R1 = sqrt(((n - m - 1) * d)^2 + a^2);
            R2 = sqrt(((n + 1 - m) * d)^2 + a^2);


            val = exp(-j * beta * R1) / R1 - 2 * exp(-j * beta * R0) / R0 + exp(-j * beta * R2) / R2;

        case 3

            psi_3_function = @(x) exp(j * beta * ((n - 1) * d + x) * cos(theta));

            val = integral(psi_3_function, 0, d);

        case 4

            psi_4_function = @(x) x .* exp(j * beta * ((n - 1) * d + x) * cos(theta));

            val = integral(psi_4_function, 0, d);

    end

end