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
function [rayTracingResults, planes] = rayTrace(antennaTrajectory, planes, satelliteTrajectories, ...
    linearizationPoint_LLA, s, startTime, samplingTime)
% RAYTRACE identifies all possible reflections between satellites and a receiver, given a list of reflector planes
% Input:
%        antennaTrajectory,         struct contaning the trajectory of the receiving antenna
%        planes,                    struct defining all planes in the urban canyon
%        satelliteTrajectories,     struct defining the satellite trajectories
%        linearizationPoint_LLA,    3x1, linearization point to convert coordinates (degrees, degrees, meters)
%        s,                         1x1, the current sample to ray trace
%        startTime,                 1x1, the starting time of the simulation (seconds)
%        samplingTime,              1x1, the sampling time of the simulation (seconds)
%
% Output:
%        rayTracingResults,         struct contaning the ray tracing results, with all identified reflection points
%        planes,                    struct contaning all the planes (this is updated since the vehicle itself is represented as several moving planes)
%
% Written by Simon Ollander
spheroid = referenceEllipsoid('WGS84');

rayTracingResults.antennaPosition_ENU = antennaTrajectory.position_ENU(:,s);
rayTracingResults.antennaVelocity_ENU = antennaTrajectory.velocity_ENU(:,s);

[antennaPosition_ECEF(1), antennaPosition_ECEF(2), antennaPosition_ECEF(3)] =...
    enu2ecef(antennaTrajectory.previousPosition_ENU(1,s),...
    antennaTrajectory.previousPosition_ENU(2,s),...
    antennaTrajectory.previousPosition_ENU(3,s),...
    linearizationPoint_LLA(1), linearizationPoint_LLA(2), ...
    linearizationPoint_LLA(3), spheroid);
rayTracingResults.antennaPosition_ECEF = antennaPosition_ECEF';

[previousAntennaPosition_ECEF(1), previousAntennaPosition_ECEF(2),...
    previousAntennaPosition_ECEF(3)] =...
    enu2ecef(antennaTrajectory.previousPosition_ENU(1,s),...
    antennaTrajectory.previousPosition_ENU(2,s),...
    antennaTrajectory.previousPosition_ENU(3,s),...
    linearizationPoint_LLA(1), linearizationPoint_LLA(2),...
    linearizationPoint_LLA(3), spheroid);
rayTracingResults.antennaVelocity_ECEF = antennaPosition_ECEF' - previousAntennaPosition_ECEF';

[antPos_LLA(1), antPos_LLA(2), antPos_LLA(3)] = enu2geodetic(...
    antennaTrajectory.position_ENU(1,s),...
    antennaTrajectory.position_ENU(2,s),...
    antennaTrajectory.position_ENU(3,s),...
    linearizationPoint_LLA(1), linearizationPoint_LLA(2), ...
    linearizationPoint_LLA(3), spheroid, 'degrees');
rayTracingResults.antennaPosition_LLA = antPos_LLA';



r = 0; % Counts reflections for each PRN
for PRN = 1:1:size(satelliteTrajectories.positions_ECEFx,1)
    
    satellitePosition_ECEF = [satelliteTrajectories.positions_ECEFx(PRN, s);...
        satelliteTrajectories.positions_ECEFy(PRN, s);...
        satelliteTrajectories.positions_ECEFz(PRN, s)];
    [satellitePosition_ENU(1), satellitePosition_ENU(2), satellitePosition_ENU(3)] = ...
        ecef2enu(satellitePosition_ECEF(1), satellitePosition_ECEF(2), satellitePosition_ECEF(3), ...
        linearizationPoint_LLA(1), linearizationPoint_LLA(2), linearizationPoint_LLA(3),...
        spheroid, 'degrees');
    
    satellitePosition_LLA = ecef2lla(satellitePosition_ECEF')';
    [satelliteAzimuth, satelliteElevation, satelliteSlantRange] = ...
        geodetic2aer(satellitePosition_LLA(1), satellitePosition_LLA(2),...
        satellitePosition_LLA(3), ...
        linearizationPoint_LLA(1), linearizationPoint_LLA(2),...
        linearizationPoint_LLA(3), spheroid, 'degrees');
    satellitePosition_AER = [satelliteAzimuth, satelliteElevation,...
        satelliteSlantRange]';
    
    satelliteVelocity_ECEF = [satelliteTrajectories.velocities_ECEFx(PRN, s); ...
        satelliteTrajectories.velocities_ECEFy(PRN, s); ...
        satelliteTrajectories.velocities_ECEFz(PRN, s)];
    
    rayTracingResults.satellitePositions_AER(1:3, PRN) = satellitePosition_AER;
    rayTracingResults.satellitePositions_ENU(1:3, PRN) = satellitePosition_ENU';
    
    rayTracingResults.satellitePositions_ECEF(1:3, PRN) = satellitePosition_ECEF;
    rayTracingResults.satelliteVelocities_ECEF(1:3, PRN) = satelliteVelocity_ECEF;
    
    lineOfSightFlag = 1;
    satelliteArrival_ENU = antennaTrajectory.position_ENU(:,s);
    for w = 1:1:length(planes)
        crossObstacle = planeLineIntersect(planes(w).normal,...
            planes(w).points(:,1), planes(w).limits, antennaTrajectory.position_ENU(:,s), satellitePosition_ENU');
        if sum(isnan(crossObstacle)) == 0
            lineOfSightFlag = 0;
            satelliteArrival_ENU = crossObstacle;
            break
        end
    end
    rayTracingResults.satelliteLineOfSightFlags(PRN) = lineOfSightFlag;
    rayTracingResults.satelliteArrivals_ENU(1:1:3, PRN) = satelliteArrival_ENU;
    
    % Compute and add reflections for this time step
    for w1 = 1:1:length(planes)
        [reflectionPoint_ENU, mirrorPoint_ENU, projectionPoint_ENU] = ...
            computeReflection(planes(w1).normal, planes(w1).points(:,1), ...
            planes(w1).limits, antennaTrajectory.position_ENU(:,s), satellitePosition_ENU');
        
        if sum(isnan(reflectionPoint_ENU)) == 0  % A reflection has been found!
            
            % Check all other walls for crossings, to see if the
            % reflection is blocked:
            reflectionBlocked = 0;
            for w2 = setdiff(1:1:length(planes),w1)
                crossObstacle_toBuilding = planeLineIntersect(planes(w2).normal,...
                    planes(w2).points(:,1), planes(w2).limits, satellitePosition_ENU', reflectionPoint_ENU);
                crossObstacle_fromBuilding = planeLineIntersect(planes(w2).normal,...
                    planes(w2).points(:,1), planes(w2).limits, reflectionPoint_ENU, ...
                    antennaTrajectory.position_ENU(:,s));
                
                % Reflection block detected, do not consider reflection!
                if  (sum(isnan(crossObstacle_toBuilding)) == 0) ...
                        || (sum(isnan(crossObstacle_fromBuilding)) == 0 )
                    reflectionBlocked = 1;
                    break
                end
                
            end
            
            if reflectionBlocked == 0 % Reflection found and it is not blocked, add to list
                r = r + 1;
                rayTracingResults.reflections(r).prn = PRN;
                
                % Save reflection points in ENU
                rayTracingResults.reflections(r).reflectionPoint_ENU = reflectionPoint_ENU;
                rayTracingResults.reflections(r).mirrorPoint_ENU = mirrorPoint_ENU;
                rayTracingResults.reflections(r).projectionPoint_ENU = projectionPoint_ENU;
                
                % Save reflection points in ECEF
                rayTracingResults.reflections(r).reflectionPoint_ECEF = enu2ecef(reflectionPoint_ENU(1),...
                    reflectionPoint_ENU(2), reflectionPoint_ENU(3),...
                    linearizationPoint_LLA(1), linearizationPoint_LLA(2),...
                    linearizationPoint_LLA(3), spheroid)';
                rayTracingResults.reflections(r).mirrorPoint_ECEF = enu2ecef(...
                    mirrorPoint_ENU(1), mirrorPoint_ENU(2), mirrorPoint_ENU(3),...
                    linearizationPoint_LLA(1), linearizationPoint_LLA(2),...
                    linearizationPoint_LLA(3), spheroid )';
                rayTracingResults.reflections(r).projPoint_ECEF = enu2ecef(projectionPoint_ENU(1), ...
                    projectionPoint_ENU(2), projectionPoint_ENU(3), linearizationPoint_LLA(1), ...
                    linearizationPoint_LLA(2), linearizationPoint_LLA(3), spheroid)';
                
                rayTracingResults.reflections(r).planeType = planes(w1).planeType;
                rayTracingResults.reflections(r).objectType = planes(w1).objectType;
                
                % Compute and check reflection angle
                vec1 = satellitePosition_ENU' - reflectionPoint_ENU; % Order here is important!
                vec2 = antennaTrajectory.position_ENU(:,s) - reflectionPoint_ENU;
                planeNorm = planes(w1).normal;
                
                inAngle = pi/2 - acos( sum(vec1.* planeNorm) / (norm(vec1)*norm(planeNorm)) );
                
                outAngle = pi/2 - acos( sum(vec2.* planeNorm) / (norm(vec2)*norm(planeNorm)) );
                
                % inAngle and outAngle relative reflector plane must be equal and between 0 and 90!
                if inAngle < 0 || inAngle > pi/2 || abs(inAngle - outAngle) > 1e-2
                    disp('Reflection logic error!')
                    disp('Check settings, something is wrong...')
                end
                rayTracingResults.reflections(r).reflectionAngle = inAngle;
                
            end
            
        end
    end
    
end
rayTracingResults.time = startTime + s*samplingTime;
rayTracingResults.noReflections = r;

if rayTracingResults.noReflections == 0
    rayTracingResults.reflections  = [];
end

end