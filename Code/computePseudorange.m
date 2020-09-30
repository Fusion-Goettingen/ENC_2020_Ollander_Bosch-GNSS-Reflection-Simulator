% Matlab software to generate an urban canyon and simulate vehicle movement to create a reflection scenario
%
% "Simulation-based analysis of multipath delay distributions in urban canyons"
% Simon Ollander, Friedrich-Wilhelm Bode and Marcus Baum
% European Navigation Conference 2020
%
% Further information:
% http://www.fusion.informatik.uni-goettingen.de
% https://github.com/Fusion-Goettingen
%
% Source code written by Simon Ollander, Bosch Car Multimedia / University of Göttingen
function [pseudoRangeError, SNR] = computePseudorange(samplingFactor, correlatorSpacing,...
    noChips_delay, noChips_SNR, phaseNoiseStandardDeviation,...
    amplitudeNoiseStandardDeviation, constellation, frequencyBands,...
    signalComponents, satellitePRN, lineOfSightFlag, noReflections, ...
    reflectionDelay, powerLossRatio, const)
% COMPUTEPSEUDORANGE takes a combination of signals as input and computes the resulting pseudorange error
% Input:
%        samplingFactor,                    1x1, describes how many times the signal should be resampled
%        correlatorSpacing,                 1x1, spacing of the correlator,between 0 and 1 (number of chips)
%        noChips_delay,                     1x1, number of chips to compute the code delay
%        noChips_SNR,                       1x1, number of chips to compute the signal-to-noise ratio
%        phaseNoiseStandardDeviation,       1x1, standard deviation of the phase noise
%        amplitudeNoiseStandardDeviation,   1x1, standard deviation of the amplitude noise
%        constellation,                     string describing the GNSS constellation ('GPS', or 'GAL')
%        frequencyBands,                    Nx1 cell string describing the frqeuency bands {'L1', 'L2'}
%        signalComponents,                  Nx1 cell string describing the signal component of each frequency ('CA')
%        satellitePRN,                      1x1, the pseudorandom number of the satellite
%        lineOfSightFlag,                   1x1, flag, 0 if the direct signal is not present, 1 if it is present
%        noReflections,                     1x1, number of received reflections
%        reflectionDelay,                   1xnoReflections, delays of all reflections
%        powerLossRatio,                    1xnoReflections, power loss ratios of all reflections
%        const,                             struct containing natural constants
%
% Output:
%        pseudoRangeError,                  Nx1 resulting pseudorange error on each frequency
%        SNR,                               Nx1 resulting SNR on each frequency
%
% Written by Simon Ollander
maxCorrelation(1) = 157542; % L1
maxCorrelation(2) = 2455200; % L2

for f = 1:1:length(frequencyBands)
    frequencyBand = frequencyBands{f};
    component = signalComponents{f};
    
    % Generate LOS time signal
    endTime = const.(constellation).(frequencyBand).(component).T_code;
    samplingFrequency = samplingFactor*const.(constellation).(frequencyBand).(component).f_carrier;
    timeVector = (0:1/samplingFrequency:endTime-1/samplingFrequency);
    code = generateAndSampleCode(constellation, frequencyBand, component,...
        satellitePRN, samplingFrequency, const);
    [signal_LOS, codeToCorrelate]...
        = generateSignal(const, constellation, frequencyBand,...
        code, timeVector, samplingFrequency, phaseNoiseStandardDeviation, ...
        amplitudeNoiseStandardDeviation, 0);
    if lineOfSightFlag == 0
        signal_LOS = zeros(1, length(timeVector));
    end
    
    % Generate multipath signal(s)
    signal_reflected_scaled = zeros(1,length(timeVector));
    for r = 1:1:noReflections
        reflectionTimeDelay = reflectionDelay(r, f)/const.c;
        [signal_reflected, ~]...
            = generateSignal(const, constellation, frequencyBand, ...
            code, timeVector, samplingFrequency,...
            phaseNoiseStandardDeviation, amplitudeNoiseStandardDeviation,...
            reflectionTimeDelay);
        signal_reflected_scaled = signal_reflected_scaled + ...
            signal_reflected .* powerLossRatio(r, f);
    end
    
    % Combine LOS and MP to complete time signal
    summedSignal = signal_LOS + signal_reflected_scaled;
    
    % Number of samples in correlation form in each direction
    noSamplesInCorrelation_SNR = noChips_SNR * samplingFrequency/const.(constellation).(frequencyBand).(component).f_chip;
    noSamplesInCorration_delay = noChips_delay * samplingFrequency/const.(constellation).(frequencyBand).(component).f_chip;
    
    % Generate correlation form (large version, for SNR computation)
    [corrForm_SNR, delays_SNR] = xcorr(...
        summedSignal.*exp(-2*1j*pi*const.(constellation).(frequencyBand).(component).f_carrier*timeVector),...
        codeToCorrelate, noSamplesInCorrelation_SNR);
    
    chipVector_SNR = delays_SNR/samplingFrequency;
    
    correlationForm_SNR_abs = abs(corrForm_SNR);
    
    % Filter out compact correlation form, for delay tracking
    centerSample_SNR = floor(length(chipVector_SNR)/2)+1;
    chipVector_delay = chipVector_SNR(centerSample_SNR - noSamplesInCorration_delay:centerSample_SNR + noSamplesInCorration_delay);
    correlationForm_delay_abs = correlationForm_SNR_abs(centerSample_SNR - noSamplesInCorration_delay:centerSample_SNR + noSamplesInCorration_delay);
    
    timeInCorrSpacing = correlatorSpacing/const.(constellation).(frequencyBand).(component).f_chip;
    
    % Method to estimate peak:
    noSamplesInCorrSpacing = timeInCorrSpacing * samplingFrequency;
    
    earlyLateDiff = NaN*ones(size(chipVector_delay));
    phaseCorrPow = NaN*ones(size(chipVector_delay));
    earlySamples = NaN*ones(size(chipVector_delay));
    phaseSamples = NaN*ones(size(chipVector_delay));
    lateSamples = NaN*ones(size(chipVector_delay));
    for s = (noSamplesInCorrSpacing + 1): 1 : (length(chipVector_delay) - noSamplesInCorrSpacing)
        
        earlySamples(s) = s - noSamplesInCorrSpacing;
        phaseSamples(s) = s;
        lateSamples(s) = s + noSamplesInCorrSpacing;
        
        % Difference between early and late:
        earlyLateDiff(s) = abs(correlationForm_delay_abs(earlySamples(s)) -...
            correlationForm_delay_abs(lateSamples(s)));
        
        % Correlation power of phase:
        phaseCorrPow(s) = correlationForm_delay_abs(phaseSamples(s));
    end
    trackingOutput = phaseCorrPow./earlyLateDiff;
    [~, bestSample] = max(trackingOutput);
    
    pseudoRangeError(f) = chipVector_delay(bestSample)*const.c;
    SNR(f) = correlationForm_delay_abs(bestSample)/maxCorrelation(f);
end


end