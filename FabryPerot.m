function result = FabryPerot(frequencies, n_air, n_sample, d, c_const, theta, beta)
sumTerm = 0;
for k = 1:2
    sumTerm = sumTerm + (FresnelReflect(n_sample, n_air, theta, beta).^2 .* ...
        waveProp_p(n_sample, frequencies, d, c_const).^2) .^ k;
end

result = 1 + sumTerm;
end

function p = waveProp_p(n, frequencies, d, c_const)
    p = exp(-1j .* n .* 2 * pi .* frequencies .* d ./ c_const);
end