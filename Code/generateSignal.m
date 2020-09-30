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
function [signal, codeToCorrelate] = generateSignal(const, constellation,...
    frequencyBand, code, timeVector, samplingFrequency,...
    phaseNoiseStandardDeviation, amplitudeNoiseStandardDeviation, reflectionTimeDelay)
% GENERATESIGNAL generates a radio-frequency signal with a given delay
% Input:
%        const,                                 struct containing natural constants
%        constellation,                         string defining the constellation 'GPS'
%        frequencyBand,                         string defining the frequency band 'L1'
%        code,                                  struct contaning the ranging code of the signal
%        timeVector,                            time vector corresponding to the signal (s)
%        samplingFrequency,                     sampling frequency of the signal (Hz)
%        phaseNoiseStandardDeviation,           standard deviation of the phase noise of the signal
%        amplitudeNoiseStandardDeviation,       standard deviation of the amplitude noise of the signal
%        reflectionTimeDelay,                   time delay of the signal with respect to the direct path (s)
%
% Output:
%        signal,                                vector contaning radio frequency signal
%        codeToCorrelate,                       struct contaning the ranging code of the signal
%
% Written by Simon Ollander

reflectionDelaySamples = round(reflectionTimeDelay*samplingFrequency);

dataBit = 1;
if strcmp(constellation, 'GPS')
    if strcmp(frequencyBand, 'L1')
        phaseDifference = 2*pi*mod(const.GPS.L1.CA.f_carrier*reflectionTimeDelay, 1);
        dataLessSignal = cos(2*pi*const.GPS.L1.CA.f_carrier*(timeVector) + phaseNoiseStandardDeviation*randn(1, length(timeVector)) + phaseDifference);
        shiftedCode = circshift(code.(constellation).(frequencyBand).CA, reflectionDelaySamples);
        signal = shiftedCode.*dataBit.*dataLessSignal + amplitudeNoiseStandardDeviation*randn(1, length(timeVector));
        codeToCorrelate = code.(constellation).(frequencyBand).CA;
    elseif strcmp(frequencyBand, 'L2')
        phaseDifference = 2*pi*mod(const.GPS.L2.CM.f_carrier*reflectionTimeDelay, 1);
        dataLessSignal = cos(2*pi*const.GPS.L2.CM.f_carrier*(timeVector) + phaseNoiseStandardDeviation*randn(1, length(timeVector)) + phaseDifference);
        shiftedCode = circshift(code.(constellation).(frequencyBand).CM, reflectionDelaySamples);
        signal = shiftedCode.*dataBit.*dataLessSignal + amplitudeNoiseStandardDeviation*randn(1, length(timeVector));
        codeToCorrelate = code.(constellation).(frequencyBand).CM;
    elseif strcmp(frequencyBand, 'L5')
        phaseDifference = 2*pi*mod(const.GPS.L5.SOL.f_carrier*reflectionTimeDelay, 1);
        dataLessSignal = cos(2*pi*const.GPS.L5.SOL.f_carrier*timeVector + phaseNoiseStandardDeviation*randn(1,length(timeVector)) + phaseDifference);
        code_L5_I = circshift(code.(constellation).(frequencyBand).I, reflectionDelaySamples);
        code_L5_Q = circshift(code.(constellation).(frequencyBand).Q, reflectionDelaySamples);
        xL5 = (code_L5_I.*dataBit + 1j*code_L5_Q);
        signal = xL5.*dataLessSignal;
        codeToCorrelate = code.(constellation).(frequencyBand).I;
    end
elseif strcmp(constellation, 'GAL')
    if strcmp(frequencyBand, 'E1')
        f_carrier = const.GAL.E1.OS.f_carrier;
        f_subcarrier = const.GAL.E1.OS.f_subcarrier;
        phaseDifference = 2*pi*mod(f_carrier*reflectionTimeDelay, 1);
        code_E1B = circshift(code.(constellation).(frequencyBand).B, reflectionDelaySamples);
        code_E1C = circshift(code.(constellation).(frequencyBand).C, reflectionDelaySamples);
        alpha = const.GAL.E1.OS.alpha;
        beta = const.GAL.E1.OS.beta;
        m = const.GAL.E1.OS.m; % 6
        n = const.GAL.E1.OS.n; % 1
        theta_MP_sc1 = 2*pi*mod(1*f_subcarrier*reflectionTimeDelay, 1);
        theta_MP_sc6 = 2*pi*mod(6*f_subcarrier*reflectionTimeDelay, 1);
        sc_E1_Ba = sign(sin(2*pi* n*f_subcarrier*timeVector + phaseNoiseStandardDeviation*randn(1,length(timeVector)) + theta_MP_sc1 ));
        sc_E1_Bb = sign(sin(2*pi* m*f_subcarrier*timeVector + phaseNoiseStandardDeviation*randn(1,length(timeVector)) + theta_MP_sc6));
        Bcomp = (code_E1B.*dataBit).*(alpha*sc_E1_Ba + beta*sc_E1_Bb);
        Ccomp = code_E1C.*1.*(alpha*sc_E1_Ba - beta*sc_E1_Bb);
        signal = ((Bcomp - Ccomp).*cos(2*pi*f_carrier*timeVector + ...
            phaseNoiseStandardDeviation*randn(1, length(timeVector)) + phaseDifference)) +...
            amplitudeNoiseStandardDeviation*randn(1, length(timeVector));
        sc = sqrt(10/11).*sign(sin(2*pi* 1*const.GAL.E1.OS.f_subcarrier*timeVector)) + sqrt(1/11).*sign(sin(2*pi*6*const.GAL.E1.OS.f_subcarrier*timeVector));
        codeToCorrelate = code.(constellation).(frequencyBand).B.*sc;
    elseif strcmp(frequencyBand, 'E5')
        phaseDifference = 2*pi*mod(const.GAL.E5.OS.f_carrier*reflectionTimeDelay, 1);
        dataLessSignal = cos(2*pi*const.GAL.E5.OS.f_carrier*timeVector + phaseNoiseStandardDeviation*randn(1,length(timeVector)) + phaseDifference);
        tbE5 = 1/const.GAL.E5.OS.f_subcarrier;
        dE5aI = dataBit;
        dE5bI = dataBit;
        eE5aI = circshift(code.(constellation).(frequencyBand).aI, reflectionDelaySamples) .* dE5aI;
        eE5aQ = circshift(code.(constellation).(frequencyBand).aQ, reflectionDelaySamples);
        eE5bI = circshift(code.(constellation).(frequencyBand).bI, reflectionDelaySamples) .* dE5bI;
        eE5bQ = circshift(code.(constellation).(frequencyBand).bQ, reflectionDelaySamples);
        
        econjE5aI = eE5aQ .* eE5bI .* eE5bQ;
        econjE5aQ = eE5aI .* eE5bI .* eE5bQ;
        econjE5bI = eE5bQ .* eE5aI .* eE5aQ;
        econjE5bQ = eE5bI .* eE5aI .* eE5aQ;
        
        AS = 1/2 * [sqrt(2)+1, 1, -1, -sqrt(2)-1, -sqrt(2)-1, -1, 1, sqrt(2)+1];
        AP = 1/2 * [-sqrt(2)+1, 1, -1, sqrt(2)-1, sqrt(2)-1, -1, 1, -sqrt(2)+1];
        
        scE5S_1per = upsampleCode(AS, tbE5*samplingFrequency/8); % /8
        scE5P_1per = upsampleCode(AP, tbE5*samplingFrequency/8); % /8
        
        delay_samples = samplingFrequency*tbE5/4;
        
        scE5S_1per_del = circshift(scE5S_1per, delay_samples);
        scE5P_1per_del = circshift(scE5P_1per, delay_samples);
        
        factor = floor(length(timeVector)/length(scE5S_1per));
        rest = length(timeVector) - factor*length(scE5S_1per);
        scE5S = [repmat(scE5S_1per, 1,factor) scE5S_1per(1:1:rest)];
        scE5S_del = [repmat(scE5S_1per_del, 1,factor) scE5S_1per_del(1:1:rest)];
        scE5P = [repmat(scE5P_1per, 1,factor) scE5P_1per(1:1:rest)];
        scE5P_del = [repmat(scE5P_1per_del, 1,factor) scE5P_1per_del(1:1:rest)];
        
        xE5 = 1/(2*sqrt(2)) * (eE5aI + 1j* eE5aQ) .* (scE5S - 1j*scE5S_del) + ...
            1/(2*sqrt(2)) * (eE5bI + 1j* eE5bQ) .* (scE5S + 1j*scE5S_del) + ...
            1/(2*sqrt(2)) * (econjE5aI + 1j* econjE5aQ) .* (scE5P - 1j*scE5P_del) + ...
            1/(2*sqrt(2)) * (econjE5bI + 1j* econjE5bQ) .* (scE5P + 1j* scE5P_del) ;
        
        signal = xE5.*dataLessSignal;
        codeToCorrelate = 1/(2*sqrt(2))*code.(constellation).(frequencyBand).aI.*(scE5S - 1j*scE5S_del);
        
    end
end

end