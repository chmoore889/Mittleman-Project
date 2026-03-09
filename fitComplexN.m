function n_fit = fitComplexN(theta_rad, measuredMag, measuredPhase, frequencies, n_guess, thickness)
%fitComplexN Fits complex N based on transfer model
%   Uses gradient descent fitting

epsilon = 0.001;
maxIterations = 1000;

prevError = Inf;
for iter = 1:maxIterations
    H_model = mittlemanTransferModel(frequencies, theta_rad, n_guess, thickness);
    
    modelMag = abs(H_model);
    modelPhase = unwrap(angle(H_model));
    modelPhase = modelPhase - modelPhase(1);
    
    magErr = measuredMag - modelMag;
    phaseErr = measuredPhase - modelPhase;
    
    totalError = sum(abs(magErr)) + sum(abs(phaseErr));
    if totalError >= prevError
        break;
    end
    prevError = totalError;
    
    %Update guess with errors
    n_guess = n_guess + epsilon .* phaseErr;%Phase error steps real component
    n_guess = n_guess - epsilon .* 1j .* magErr;%Mag error steps imag component
end

if iter == maxIterations
    warning("Limited by iterations")
end

n_fit = n_guess;
end