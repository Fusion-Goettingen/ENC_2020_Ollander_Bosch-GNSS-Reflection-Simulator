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
function resampledCodes = generateAndSampleCode(constellation, frequencyBand, component,...
    satellitePRN, samplingFrequency, const)
% GENERATEANDSAMPLECODE generates code for a satellite and resamples it to match the sampling frequency
% Input:
%        constellation,     string with the GNSS constellation, ('GPS')
%        frequencyBand,     string with frequency band, ('L1')
%        component,         string with the component of the signal ('CA')
%        satellitePRN,      1x1, the pseudorandom number of the satellite
%        samplingFrequency, 1x1, the sampling frequency of the signal (Hz)
%        const,             struct containing natural constants
%
% Output:
%        resampledCodes,    struct containing the satellite codes, resampled to match the sampling frequency
%
% Written by Simon Ollander

if strcmp(constellation, 'GPS')
    if strcmp(frequencyBand, 'L1')
        load('codes_L1CA.mat')
        codes.(constellation).(frequencyBand).CA = codes_L1CA(:, satellitePRN)';
    elseif strcmp(frequencyBand, 'L2')
        if strcmp(component, 'CL')
            load('codes_L2CL.mat')
            codes.(constellation).(frequencyBand).CL = codes_L2CL(:, satellitePRN)';
        elseif strcmp(component, 'CM')
            load('codes_L2CM.mat')
            codes.(constellation).(frequencyBand).CM = codes_L2CM(:, satellitePRN)';
        end
    elseif strcmp(frequencyBand, 'L5')
        load('codes_L5I.mat')
        load('codes_L5Q.mat')
        codes.(constellation).(frequencyBand).I = codes_L5I(:, satellitePRN)';
        codes.(constellation).(frequencyBand).Q = codes_L5Q(:, satellitePRN)';
    end
elseif strcmp(constellation, 'GAL')
    if strcmp(frequencyBand, 'E1')
        load('codes_E1B.mat')
        load('codes_E1C.mat')
        codes.(constellation).(frequencyBand).B = codes_E1B(:, satellitePRN)';
        codes.(constellation).(frequencyBand).C = codes_E1C(:, satellitePRN)';
    elseif strcmp(frequencyBand, 'E5')
        load('codes_E5aI.mat')
        load('codes_E5aQ.mat')
        load('codes_E5bI.mat')
        load('codes_E5bQ.mat')
        codes.(constellation).(frequencyBand).aI = codes_E5aI(:, satellitePRN)';
        codes.(constellation).(frequencyBand).aQ = codes_E5aQ(:, satellitePRN)';
        codes.(constellation).(frequencyBand).bI = codes_E5bI(:, satellitePRN)';
        codes.(constellation).(frequencyBand).bQ = codes_E5bQ(:, satellitePRN)';
    end
end

allFields = fields(codes.(constellation).(frequencyBand));
for f = 1:1:length(allFields)
    code = codes.(constellation).(frequencyBand).(allFields{f});
    if strcmp(constellation, 'GAL')
        setLength = const.(constellation).(frequencyBand).OS.T_code*samplingFrequency;
    elseif strcmp(constellation, 'GPS') && strcmp(frequencyBand, 'L5')
        const.(constellation).(frequencyBand).SOL.T_code;
        setLength = const.(constellation).(frequencyBand).SOL.T_code*samplingFrequency;
    else
        setLength = const.(constellation).(frequencyBand).(allFields{f}).T_code*samplingFrequency;
    end
    
    % If needed, resample code
    resamplingFactor = setLength/size(code, 2);
    
    % code_upsampled = upsampleCode(currCode, factor);
    if resamplingFactor ~= 1
        codeLength = length(code);
        % fractional upsampling, zero order hold
        index = 0;
        for count = 1/resamplingFactor:1/resamplingFactor:codeLength
            index = index + 1;
            if ceil(count) > codeLength
                resampledCode_temp(:, index) = code(:, codeLength);
            else
                resampledCode_temp(:, index) = code(:, ceil(count));
            end
        end
        resampledCode = resampledCode_temp;
    else
        resampledCode = code;
    end
    
    resampledCodes.(constellation).(frequencyBand).(allFields{f}) = resampledCode;
end

end