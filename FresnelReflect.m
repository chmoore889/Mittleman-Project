function reflect = FresnelReflect(n_a, n_b, theta, beta)
numerator = n_b .* cos(theta) - n_a .* cos(beta);
denominator = n_a .* cos(beta) + n_b .* cos(theta);
reflect = numerator ./ denominator;
end