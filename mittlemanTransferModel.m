function H_omega = mittlemanTransferModel(frequencies, theta_rad, n_sample, thickness)
n_air = 1.00027;
c = 299792458; %m/s

beta = asin(n_air .* sin(theta_rad) ./ n_sample);
d = thickness ./ cos(beta);

%assert(theta_rad >= beta);
m = d .* cos(theta_rad - beta);

firstTermNum = 4 .* n_air .* n_sample .* cos(theta_rad) .* cos(beta);
firstTermDenom = (n_air .* cos(beta) + n_sample .* cos(theta_rad)).^2;
firstTerm = firstTermNum ./ firstTermDenom;

exponentialNum = -1j * (d .* n_sample - m .* n_air) .* 2 * pi .* frequencies;
exponentialTerm = exp(exponentialNum ./ c);

H_omega = firstTerm .* exponentialTerm .* FabryPerot(frequencies, n_air, n_sample, d, c, theta_rad, beta);
end