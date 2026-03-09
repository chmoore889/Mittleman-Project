function [frequency, transformResult] = freqeuncyTransform(samplePeriodPs, THz)
    L = length(THz);
    L = 2^nextpow2(L);

    Y = fft(THz, L);
    P2 = Y/L;
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    transformResult = P1;

    frequency = (0:L/2) / (L * samplePeriodPs * 1e-12);
end