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
function [xPositions, yPositions, zPositions, xVelocities, yVelocities, zVelocities, prnLegend]...
    = generateSatelliteTrajectories(startTime, samplingTime, duration, constellations, const)
% GENERATESATELLITETRAJECTORIES generates satellite positions and velocities for a given time interval
% Input:
%        startTime,         1x1, start time of the interval (seconds)
%        samplingTime,      1x1, sampling time of the satellite trajectories (seconds)
%        duration,          1x1, duration of the satellite trajectories (seconds)
%        constellations,    cell string defining the GNSS constellations to be used {'GPS', 'GAL'}
%        const,             struct containing natural constants
%
% Output:
%        xPositions,        ECEF x positions of the satellite trajectories (meters)
%        yPositions,        ECEF y positions of the satellite trajectories (meters)
%        zPositions,        ECEF z positions of the satellite trajectories (meters)
%        xVelocities,       ECEF x velocities of the satellite trajectories (meters/second)
%        yVelocities,       ECEF y velocities of the satellite trajectories (meters/second)
%        zVelocities,       ECEF z velocities of the satellite trajectories (meters/second)
%        prnLegend,         list of the order of the satellites, to keep track of the constellation
%
% Written by Simon Ollander
%Based on: https://de.mathworks.com/matlabcentral/fileexchange/65845-satellite-constellation

earthGraviation = 3.986004418e+14; % Earth gravitational parameter (m^3/sec^2)

xPositions = [];
yPositions = [];
zPositions = [];

xVelocities = [];
yVelocities = [];
zVelocities = [];

prnLegend = [];

for c = 1:1:length(constellations)
    
    constellation = constellations{c};
    
    semiMajorAxis = const.orbits.(constellation).semiMajorAxis;
    meanAnomaly = const.orbits.(constellation).meanAnomaly;
    perigeeArgument = const.orbits.(constellation).perigeeArgument;
    orbitalInclination = const.orbits.(constellation).orbitalInclination;
    rightAscensionAscendingNode = const.orbits.(constellation).rightAscensionAscendingNode;
    eccentricity = const.orbits.(constellation).eccentricity;
    prnOrder = const.orbits.(constellation).prnOrder;
    
    
    numberOfSatellites = length(semiMajorAxis);
    clear currXs currYs currZs
    for s = 1:1:numberOfSatellites
        globalPRN = s;
        if strcmp(constellation, 'GAL')
            globalPRN = globalPRN + 31;
        end
        prnLegend = [prnLegend, globalPRN];
        
        timeVector = startTime:samplingTime:(startTime + duration);
        
        xPositions_temp = zeros(1, length(timeVector));
        yPositions_temp = zeros(1, length(timeVector));
        zPositions_temp = zeros(1, length(timeVector));
        
        for step = 1:1:length(timeVector)
            
            Mi = meanAnomaly(1, s) + sqrt(earthGraviation/semiMajorAxis(1, s)^3)*timeVector(step);
            
            E1 = Mi;
            E2 = E1-((E1-eccentricity(1, s)*(sin(E1))-Mi)/(1-eccentricity(1, s)*(cos(E1))));
            
            tolerance = 0.001;
            error = 100;
            
            while error > tolerance
                E1 = E2;
                E2 = E1-((E1-eccentricity(1, s)*(sin(E1))-Mi)/(1-eccentricity(1, s)*(cos(E1))));
                error = abs(E2-E1);
            end
            radius = semiMajorAxis(1, s)*(1-eccentricity(1, s)*cos(E2));
            trueAnomaly = 2*atan(sqrt((1-eccentricity(1, s))/(1+eccentricity(1, s)))*tan(E2/2));
            u = perigeeArgument(1, s) + trueAnomaly;
            
            X = radius*cos(u);
            Y = radius*sin(u);
            
            xPositions_temp(:, step) = X*cos(rightAscensionAscendingNode(1, s))-Y*cos(orbitalInclination(1, s))*sin(rightAscensionAscendingNode(1, s)); % ECEF x-coordinate SAT (meters)
            yPositions_temp(:, step) = X*sin(rightAscensionAscendingNode(1, s))+Y*cos(orbitalInclination(1, s))*cos(rightAscensionAscendingNode(1, s)); % ECEF y-coordinate SAT (meters)
            zPositions_temp(:, step) = Y*sin(orbitalInclination(1, s)); % ECEF z-coordinate SAT (meters)
        end
        satelliteXPosition(s, :) = [xPositions_temp];
        satelliteYPosition(s, :)= [yPositions_temp];
        satelliteZPosition(s, :)= [zPositions_temp];
        
    end
    
    % Satellite velocities
    currXv = diff(satelliteXPosition, 1, 2)/samplingTime;
    currYv = diff(satelliteYPosition, 1, 2)/samplingTime;
    currZv = diff(satelliteZPosition, 1, 2)/samplingTime;
    
    % Remove first sample to get same amount of samples as in velocities
    satelliteXPosition(:, 1) = [];
    satelliteYPosition(:, 1) = [];
    satelliteZPosition(:, 1) = [];
    
    xPositions = [xPositions; satelliteXPosition];
    yPositions = [yPositions; satelliteYPosition];
    zPositions = [zPositions; satelliteZPosition];
    
    xVelocities = [xVelocities; currXv];
    yVelocities = [yVelocities; currYv];
    zVelocities = [zVelocities; currZv];
end

end