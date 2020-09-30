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
function D = computePhysics(rejectionLHCP, carrierFrequencies, D, const)
% COMPUTEPHYSICS description
% Input:
%        rejectionLHCP,         left hand circular polarization rejection of the antenna (dB)
%        carrierFrequencies,    1xN, vector containing N carrier frequencies (Hz)
%        D,                     struct containing all reflections
%        const,                 struct containing natural constants
%
% Output:
%        D,                     struct containing all reflections, including delay, power loss, Doppler shift, and carrier phase difference
%
% Written by Simon Ollander

noFreqencies = length(carrierFrequencies);
spheroid = referenceEllipsoid('WGS84');

for d = 1:1:length(D)
    for c = 1:1:size(D(d).data.results,1)
        for s = 1:1:D(d).data.noSamples
            results = D(d).data.results(c, s);
            for r = 1:1:results.noReflections
                reflection = results.reflections(r);
                PRN = results.reflections(r).prn;
                
                % Compute delay of reflection
                reflectionDelay = norm(results.antennaPosition_ENU - reflection.reflectionPoint_ENU)...
                    + norm(results.satellitePositions_ENU(:, PRN) - reflection.reflectionPoint_ENU)...
                    - norm(results.antennaPosition_ENU - results.satellitePositions_ENU(:, PRN));
                D(d).data.results(c, s).reflections(r).reflectionDelay(1:noFreqencies) = ...
                    reflectionDelay;
                
                % Compute 3D Doppler of direct signal
                antennaPosition_ECEF = lla2ecef(results.antennaPosition_LLA')';
                
                [results.antennaVelocity_ECEF(1), results.antennaVelocity_ECEF(2), results.antennaVelocity_ECEF(3)]...
                    = enu2ecefv(results.antennaVelocity_ENU(1), results.antennaVelocity_ENU(2),...
                    results.antennaVelocity_ENU(3),...
                    D.data.initialAntennaPosition_LLA(1), D.data.initialAntennaPosition_LLA(2));
                frequencyLineOfSight = computeDoppler(antennaPosition_ECEF, results.antennaVelocity_ECEF',...
                    results.satellitePositions_ECEF(:, PRN), ...
                    results.satelliteVelocities_ECEF(:, PRN), ...
                    carrierFrequencies, const);
                
                [reflectionPoint_ECEF(1), reflectionPoint_ECEF(2), reflectionPoint_ECEF(3)] =...
                    enu2ecef(reflection.reflectionPoint_ENU(1), reflection.reflectionPoint_ENU(2),...
                    reflection.reflectionPoint_ENU(3),...
                    D.data.initialAntennaPosition_LLA(1), D.data.initialAntennaPosition_LLA(2),...
                    D.data.initialAntennaPosition_LLA(3), spheroid);
                frequencyToBuilding = computeDoppler(reflectionPoint_ECEF', [0, 0, 0]',...
                    results.satellitePositions_ECEF(:,PRN), ...
                    results.satelliteVelocities_ECEF(:,PRN), carrierFrequencies, const);
                
                frequencyReflection = computeDoppler(antennaPosition_ECEF, results.antennaVelocity_ECEF',...
                    reflectionPoint_ECEF', [0, 0, 0]', frequencyToBuilding, const);
                dopplerDifference = frequencyReflection - frequencyLineOfSight;
                
                D(d).data.results(c, s).reflections(r).dopplerDifference(1:noFreqencies)...
                    = dopplerDifference;
                
                % Decide material
                material = getMaterial(results.reflections(r).objectType);
                D(d).data.results(c, s).reflections(r).material = material;
                
                reflectionCoefficient = computeReflectionCoefficient(const, material, frequencyToBuilding, ...
                    reflection.reflectionAngle, rejectionLHCP, results.satellitePositions_AER(2, PRN));
                
                % Power loss
                powerLossRatio = abs(reflectionCoefficient.Effective);
                powerLoss = -20*log10(powerLossRatio);
                
                D(d).data.results(c, s).reflections(r).powerLoss(1:noFreqencies) = powerLoss; % decibel
                
                % Phase difference
                phaseDifference = angle(reflectionCoefficient.Effective);
                D(d).data.results(c, s).reflections(r).phaseDifference(1:noFreqencies) = phaseDifference; % multi-freq!
                
                
            end
        end
    end
end

end