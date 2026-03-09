function transmit = FresnelTransmit(n_a, n_b, theta, beta)
numerator = 2 .* n_a .* cos(theta);
denominator = n_a .* cos(beta) + n_b .* cos(theta);
transmit = numerator ./ denominator;
end