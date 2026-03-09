clear;
close all;

%% Data import
referenceFile = "samples/ref_Samp2.txt";
sampleFile = "samples/samp_Samp2.txt";

[refTimePs, refTHz] = importTHzFile(referenceFile);
[sampleTimePs, sampleTHz] = importTHzFile(sampleFile);
assert(length(refTimePs) == length(sampleTimePs));

%Plot raw time traces
figure;
plot(refTimePs, refTHz, 'DisplayName', 'Reference');
hold on;
plot(sampleTimePs, sampleTHz, 'DisplayName', 'Sample');
title("Time Series");
xlabel("Time (ps)");
ylabel("E");
legend;

timePeriodPs = mean(diff(refTimePs));
assert(abs(timePeriodPs - mean(diff(sampleTimePs))) < 1e-6);

%Find best delay for pulse overlap for pulse time delay
startOffsetPs = sampleTimePs(1) - refTimePs(1);
[corr, lags] = xcorr(sampleTHz, refTHz, 'unbiased');
[~, idx] = max(abs(corr));
pulseDelayPs = lags(idx) * timePeriodPs;
pulseDelayPs = pulseDelayPs + startOffsetPs;%Compensate if start times are not same

disp("Loaded time traces");
clearvars -except timePeriodPs refTHz sampleTHz pulseDelayPs;
%% Fourier Transform
frequencyRange = [0.25e12 1.8e12];
[refFrequencies, refTransform] = freqeuncyTransform(timePeriodPs, refTHz);
[sampleFrequencies, sampleTransform] = freqeuncyTransform(timePeriodPs, sampleTHz);

%Plot frequency spectrum magnitudes
figure;
plot(refFrequencies, abs(refTransform), 'DisplayName', 'Reference');
hold on;
plot(sampleFrequencies, abs(sampleTransform), 'DisplayName', 'Sample');
title("Reference FFT Magnitude");
xlabel("Frequency (Hz)");
ylabel("|FFT|");
xlim([0 2e12]);
xregion(frequencyRange, 'DisplayName', 'Deconvolution range');
legend;

assert(all(refFrequencies == sampleFrequencies));
frequencies = refFrequencies;

disp("Transformed to frequency domain");
clearvars -except frequencies refTransform sampleTransform timePeriodPs pulseDelayPs frequencyRange;
%% Deconvolution
indexSelection = frequencies <= frequencyRange(2) & frequencies >= frequencyRange(1);
% assert(all(abs(sampleTransform(indexSelection)) < abs(refTransform(indexSelection))));

H_omega = sampleTransform ./ refTransform;
H_omega(~indexSelection) = 0;

%Plot inverse transform of deconvolution
Y2 = [H_omega(1); H_omega(2:end)/2; flipud(conj(H_omega(2:end)))/2];
deconvolveTHz = ifft(Y2);

figure;
plot((0:length(deconvolveTHz)-1) * timePeriodPs, deconvolveTHz);
title("Deconvolution");
xlabel("Time (ps)");
ylabel("E");

%Only keep frequencies we will use
frequencies = frequencies(indexSelection);
H_omega = H_omega(indexSelection);
frequencies = frequencies(:);
H_omega = H_omega(:);

disp("Deconvolved reference from sample");
clearvars -except frequencies H_omega pulseDelayPs;
%% Fitting
n1 = 1.2; %Low limit on index
n2 = 8; %High limit on index
n_air = 1.00027;
c = 299792458; %m/s

theta_deg = 0; % Sample angle in degrees
theta_rad = deg2rad(theta_deg);

l_upper = pulseDelayPs * 1e-12 * c / (n1 - n_air);
l_lower = pulseDelayPs * 1e-12 * c / (n2 - n_air);

%Get magnitude, unwrap phase, set initial phase to 0
measuredMag = abs(H_omega);
measuredPhase = unwrap(angle(H_omega));
measuredPhase = measuredPhase - measuredPhase(1);

%Start fitting passes
lastPass = 3;
for passIndex = 1:lastPass
    thicknesses = linspace(l_lower, l_upper, 100);
    TV_or_2 = zeros(length(thicknesses), 1);%Used as TV for initial passes TV2 for pass 3
    
    for i = 1:length(thicknesses)        
        %Initial guesses for n using assumption of flat frequency response
        n_guess = (pulseDelayPs * 1e-12 * c / thicknesses(i)) + n_air;
        n_guess = complex(n_guess * ones(size(frequencies)));
        
        %Gradient descent fit for n
        n_fit = fitComplexN(theta_rad, measuredMag, measuredPhase, frequencies, n_guess, thicknesses(i));
        
        %Calculate TV or TV2
        D = diff(n_fit);
        D = abs(real(D)) + abs(imag(D));

        if passIndex == lastPass
            %TV2 calculation
            if i > 1
                TV_or_2(i) = sum(abs(diff(D)));
            end
        else
            %TV calculation
            TV_or_2(i) = sum(D);
        end
    end

    if passIndex == lastPass
        %When in the 3rd pass and using TV2, first element should be 0
        assert(TV_or_2(1) == 0);
        TV_or_2(1) = [];

        %Since differences of D are used for TV2, interpolate adjacent
        %thicknesses for better estimation of thickness
        thicknesses = thicknesses(1:end-1) + diff(thicknesses)/2;
        assert(length(thicknesses) == length(TV_or_2));
    end
    
    %Find local minima and pick the deepest
    [~, localMins, ~, prominences] = findpeaks(-TV_or_2); %findpeaks on negated to get minima
    if isempty(localMins)
        error("No local minima found in pass %d", passIndex);
    end
    [~, deepest] = max(prominences);%Use most prominent local minima
    %[~, deepest] = min(TV_or_2(localMins));
    bestThickness = thicknesses(localMins(deepest));
    thicknessStep = mean(diff(thicknesses));

    %Set new thickness bounds based on best thickness
    l_lower = bestThickness - thicknessStep * 10;
    l_upper = bestThickness + thicknessStep * 10;

    %Plot pass results
    figure;
    plot(thicknesses, TV_or_2);
    hold on;
    xregion(l_lower, l_upper);
    plotTitle = sprintf("Pass %d", passIndex);
    title(plotTitle);
    xlabel("Thickness (m)");
    if passIndex == lastPass
        ylabel("TV2");
    else
        ylabel("TV");
    end
    
    fprintf("Finished pass #%d\n", passIndex);
end

%Get final n based on thickness
n_guess = (pulseDelayPs * 1e-12 * c / bestThickness) + n_air;
n_guess = complex(n_guess * ones(size(frequencies)));
n_fit = fitComplexN(theta_rad, measuredMag, measuredPhase, frequencies, n_guess, bestThickness);

%Plot final results
thickness_mm = bestThickness * 1e3;
thicknessStep_mm = thicknessStep * 1e3;
thicknessTitle = sprintf("Thickness (mm): %f ± %f", thickness_mm, thicknessStep_mm);

figure;
frequenciesTHz = frequencies * 1e-12;
subplot(2, 1, 1);
plot(frequenciesTHz, real(n_fit));
title("n(\omega)");
xlabel("\omega(THz)");
xlim([min(frequenciesTHz) max(frequenciesTHz)]);
subtitle(thicknessTitle);

subplot(2, 1, 2);
plot(frequenciesTHz, -imag(n_fit))
title("\kappa(\omega)");
xlabel("\omega(THz)");
xlim([min(frequenciesTHz) max(frequenciesTHz)]);
subtitle(thicknessTitle);

%Figure 3 replication
figure;
l_upper = pulseDelayPs * 1e-12 * c / (n1 - n_air);
l_lower = pulseDelayPs * 1e-12 * c / (n2 - n_air);
thicknessesToShow = linspace(l_lower, l_upper, 6);
thicknessesToShow(end + 1) = bestThickness;
thicknessesToShow = sort(thicknessesToShow);

for i = 1:length(thicknessesToShow)
    n_guess = (pulseDelayPs * 1e-12 * c / thicknessesToShow(i)) + n_air;
    n_guess = complex(n_guess * ones(size(frequencies)));
    nToShow = fitComplexN(theta_rad, measuredMag, measuredPhase, frequencies, n_guess, thicknessesToShow(i));

    subplot(2, 1, 1);
    plot(frequencies, real(nToShow), 'DisplayName', sprintf("Thickness %f mm", thicknessesToShow(i) * 1e3));
    hold on;

    subplot(2, 1, 2);
    plot(frequencies, -imag(nToShow), 'DisplayName', sprintf("Thickness %f mm", thicknessesToShow(i) * 1e3));
    hold on;
end
title("\kappa");
legend;
subplot(2, 1, 1);
title("n")
legend;