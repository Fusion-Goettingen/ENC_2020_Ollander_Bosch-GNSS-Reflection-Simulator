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
%%
% Part 1: Define and generate an urban canyon
init;

tempDataFolder = '../Data_temp'; % Define where temporary data should be stored

geometry.block.size.x = 250; % Size of building block in x direction (meters)
geometry.block.size.y = 250; % Size of building block in x direction (meters)
geometry.block.noBlocks.x = 3; % Number of blocks in x direction
geometry.block.noBlocks.y = 3; % Number of blocks in y direction

geometry.build.meanHeight = 30; % Mean height of each building along the road (meters)
geometry.build.deviationHeight = 5; % Standard deviation of height of each building along the road (meters)

geometry.build.meanWidth = 25; % Mean width of each building along the road (meters)
geometry.build.deviationWidth = 0; % Standard deviation of building width (meters)

geometry.build.depth = 65; % Depth of each building, distance between road and inner court yard (meters)
geometry.build.maxDepthOffset = 0; % Maximal distance between road and building edge, uniform from 0 to max value (meters)

geometry.road.width = 30; % Width of the road (meters)

geometry.skipCentralBlock = 0; % Set variable to 0 to create a continuous urban canyon, 1 to skip the central building block to create an open square scenario
geometry.centerShift = [0, 0, 0]'; % Shifts the centerpoint in x, y, z directions (meters)

disp('Generating urban canyon...')
[planes, blockCorners] = generateCity(geometry);
noGeoPlanes = length(planes);
disp('Urban canyon generated!')
%%
% Part 2: Visualize the generated urban canyon

font = 'latex' % Interpreter for the text in the plots
fontSize = 30 % Font size of the text in the plots
plotLineWidth = 2

zoomPoint = [450, 450, 0]' % Select x, y, z point to zoom in around (meters)
zoomSize = 1000; % Select how large volume around the zoom point that should be displayed (meters)

set(0, 'DefaultTextInterpreter', font)
set(0, 'DefaultLegendInterpreter', font)
set(0, 'DefaultAxesTickLabelInterpreter', font)
set(0, 'defaultTextFontSize', fontSize)
set(0, 'defaultAxesFontSize', fontSize)
set(0, 'defaultlinelinewidth', plotLineWidth);


% Define the colors of the objects
planeColors.car = 0.8*ones(1,3);
planeColors.ground = 0.5*ones(1,3);
planeColors.build = 0.3*ones(1,3);

close all
figure('units', 'normalized', 'outerposition', [0 0 1 1])
box on
hold on

% Plot all planes
for w = 1:1:length(planes)
    plotRectangle(planes(w).points(:,1), planes(w).points(:,2),...
        planes(w).points(:,3), planes(w).points(:,4),...
        planeColors.(planes(w).objectType))
end

view(-10, 43)
xlabel('East (m)')
ylabel('North (m)')
zlabel('Up (m)')

xlim([zoomPoint(1)-zoomSize/2, zoomPoint(1)+zoomSize/2])
ylim([zoomPoint(2)-zoomSize/2, zoomPoint(2)+zoomSize/2])
zlim([0 zoomSize])
zlim([0 50])
%%
% Part 3: Define vehicle size and trajectory
constellationForwarding = 720; % Time jump of each constellation fast-forward (minutes)
samplingTime = 4; % Sampling time of the scenario (seconds)

car.height = 1.5; % Car height (meters)
car.length = 2; % Car length (meters)
car.width = 2; % Car width (meters)
car.distanceAntennaRoof = 0.01; % Distance between highest point of antenna and car roof (meters)

initialAntennaLatitude = 40; % Initial latitude of the antenna (degrees)
initialAntennaLongitude = -73; % Initial longitude of the antenna (degrees)

initialAntennaVelocity_ENU = [0, 5, 0]'; % Initial velocity vector of the antenna (meters per second, east , north, up)

trajectoryRoadCentrality = 1/2; % Define how far to the right the car should drive (1 = at the right edge, 0 = at the left edge)

initialAntennaEastPosition = blockCorners.x(2) - geometry.road.width*trajectoryRoadCentrality;
initialAntennaNorthPosition = geometry.block.size.y + geometry.road.width*1.5; 

antennaAltitude = car.height + car.distanceAntennaRoof;
initialAntennaPosition_ENU = [initialAntennaEastPosition, initialAntennaNorthPosition, antennaAltitude]';
initialAntennaPosition_LLA = [initialAntennaLatitude, initialAntennaLongitude, antennaAltitude]';

linearizationPoint_LLA = initialAntennaPosition_LLA;

constellationStartTimes = 1:(100*(constellationForwarding*60)/86400):50;

distancePerTurn = geometry.block.size.y + geometry.road.width;
sampelsPerHeading = distancePerTurn/(norm(initialAntennaVelocity_ENU)*samplingTime);

repeatedHeadings = repmat(0:90:270, sampelsPerHeading, 1);
headings = repeatedHeadings(:)';
simulationDuration = (distancePerTurn*4)/norm(initialAntennaVelocity_ENU);

samplingFrequency = 1/samplingTime;
noSamples = floor(simulationDuration*samplingFrequency);
noConstellationRepetitions = length(constellationStartTimes);
%noSteps = noSamples * noConstellationRepetitions;

previousAntennaPosition_ENU = initialAntennaPosition_ENU;
clear antennaTrajectory
for s = 1:1:noSamples
    antennaHeading_radians = headings(s)*pi/180;
    antennaVelocity_ENU = inv([cos(antennaHeading_radians), -1*sin(antennaHeading_radians), 0;
        sin(antennaHeading_radians), cos(antennaHeading_radians), 0;
        0, 0, 1])*initialAntennaVelocity_ENU;
    
    antennaPosition_ENU = previousAntennaPosition_ENU + samplingTime*antennaVelocity_ENU;
    
    antennaTrajectory.position_ENU(1:1:3, s) = antennaPosition_ENU;
    antennaTrajectory.previousPosition_ENU(1:1:3, s) = previousAntennaPosition_ENU;
    antennaTrajectory.velocity_ENU(1:1:3, s) =...
        (antennaPosition_ENU - previousAntennaPosition_ENU)/samplingTime;
    
    previousAntennaPosition_ENU = antennaPosition_ENU;
end
hold on 
plot(antennaTrajectory.position_ENU(1,:), antennaTrajectory.position_ENU(2,:))
%%
% Part 4: Ray tracing (linear algebra computation of reflections)
% This is the most time-consuming part

constellationsToUse = {'GPS'}; % List of constellations to simulate. GPS, Galileo, or both

disp('Starting computation of reflections...')
clear results satelliteTrajectories
computationTimes = [];
for c = 1:1:noConstellationRepetitions
    startTime = 86400/samplingFrequency * constellationStartTimes(c) / 100
    
    [satelliteTrajectories.positions_ECEFx, satelliteTrajectories.positions_ECEFy,...
    satelliteTrajectories.positions_ECEFz, satelliteTrajectories.velocities_ECEFx,...
    satelliteTrajectories.velocities_ECEFy, satelliteTrajectories.velocities_ECEFz] =...
        generateSatelliteTrajectories(startTime, samplingTime, simulationDuration, constellationsToUse, const);
    
    for s = 1:1:noSamples
        tic
        
        % This works if car is symmetric (otherwise car rotation must be
        % added)
        carUnderSide = antennaTrajectory.position_ENU(:,s) - 1*[car.width/2, car.length/2, antennaAltitude]';
        carUpperSide = antennaTrajectory.position_ENU(:,s)+ [car.width/2, car.length/2, car.height-antennaAltitude]';
        carPlane = generateBlock(carUnderSide, carUpperSide, 'car'); % Last slot is always the car
        planes(noGeoPlanes+1:end) = [];
        planes = [planes, carPlane];
        
        [results(c, s), planes] = rayTrace(antennaTrajectory, planes, satelliteTrajectories, ...
            linearizationPoint_LLA, s, startTime, samplingTime);
        
        computationTimes = [computationTimes, toc];
        if mod(s,10) == 0 % Output time estimation every 10th sample
            remainingSteps = noSamples*(noConstellationRepetitions - c) + noSamples - s;
            remainingMinutes = mean(computationTimes) * remainingSteps / 60;
            disp([num2str(remainingMinutes,2) ' minutes remaining'])
        end    
    end  
end

disp('All reflections generated!')
save([tempDataFolder '/' 'rayTracingSimulation' num2str(geometry.build.meanHeight) '.mat'])
disp('Results saved!')
%%
% Part 5: When part 4 is done for all wanted building configurations, save
% all the results into one aggregated struct

init
tempDataFolder = '../Data_temp';
meanBuildingHeights = [30] % List of all raytraced mean building heights (meter)

numberOfSimulations = length(meanBuildingHeights)
disp('Loading simulation results...')
clear D
for d = 1:1:numberOfSimulations
    disp(['Data ' num2str(d/numberOfSimulations*100) '% loaded'])
    
    D(d).data = load([tempDataFolder '/' 'rayTracingSimulation' num2str(meanBuildingHeights(d)) '.mat']);
end

save([tempDataFolder '/aggregatedResults.mat'])
disp('Aggregated simulation results saved!')
%%
% Part 6: Load the aggregated struct
% Simulation can be restarted from here, so the previous steps don't have
% to be resimulated

init
tempDataFolder = '../Data_temp';
load([tempDataFolder '/aggregatedResults.mat'])
disp('Data loaded!')
%%
% Part 7: Compute reflections (physics part)
rejectionLHCP = 3; % Antenna rejection of left hand circular polarized signals (dB)
carrierFrequencies = [const.GPS.L1.CA.f_carrier, const.GPS.L2.CL.f_carrier,...
    const.GPS.L5.SOL.f_carrier]; % list of carrier frequencies to simulate (Hz)
disp('Computing physics...')
D = computePhysics(rejectionLHCP, carrierFrequencies, D, const);
disp('Physics computed!')
%%
% Part 8: Aggregate all reflections to compute their distributions
disp('Computing distributions...')
D = computeDistribution(D);
disp('Distributions computed!')
%%
% Part 9: Visualize urban canyon and raytracing results

font = 'latex' % Interpreter for the text in the plots
fontSize = 30 % Font size of the text in the plots
plotLineWidth = 2

zoomPoint = [450, 450, 0]' % Select x, y, z point to zoom in around (meters)
zoomSize = 1000; % Select how large volume around the zoom point that should be displayed (meters)

simulationToPlot = 1
constellationRepetitionToPlot = 1
sampleToPlot = 50

set(0, 'DefaultTextInterpreter', font)
set(0, 'DefaultLegendInterpreter', font)
set(0, 'DefaultAxesTickLabelInterpreter', font)
set(0, 'defaultTextFontSize', fontSize)
set(0, 'defaultAxesFontSize', fontSize)
set(0, 'defaultlinelinewidth', plotLineWidth);

% Attribute one color to visualuze each received satellite
possibleColors = {'r', 'g', 'b', 'm', 'c', 'k', 'y', [255,192,203]/255,...
    [0.5 0.5 0.5], [255, 165, 0]/255, [0.1 0.1 0.1]};

% Define the colors of the objects
planeColors.car = 0.8*ones(1,3);
planeColors.ground = 0.5*ones(1,3);
planeColors.build = 0.3*ones(1,3);

close all
figure('units', 'normalized', 'outerposition', [0 0 1 1])
box on
hold on

data = D(simulationToPlot).data;
results = data.results(constellationRepetitionToPlot, sampleToPlot);

% Plot all planes
for w = 1:1:length(data.planes)-5
    plotRectangle(data.planes(w).points(:,1), data.planes(w).points(:,2),...
        data.planes(w).points(:,3), data.planes(w).points(:,4),...
        planeColors.(data.planes(w).objectType))
end

view(-10, 43)
xlabel('East (m)')
ylabel('North (m)')
zlabel('Up (m)')

xlim([zoomPoint(1)-zoomSize/2, zoomPoint(1)+zoomSize/2])
ylim([zoomPoint(2)-zoomSize/2, zoomPoint(2)+zoomSize/2])
zlim([0 zoomSize])
zlim([0 50])

% Plot reflection(s)
disp('Plotting signal paths...')
numberOfReceivedSatelltes = length(D(simulationToPlot).receivedReceptionModes(constellationRepetitionToPlot).list);
colorsToUse = possibleColors(1:1:numberOfReceivedSatelltes);

hold on
scatter3(results.antennaPosition_ENU(1), results.antennaPosition_ENU(2),...
    results.antennaPosition_ENU(3), 100, 'k', 'Marker', 'o') 

for s = 1:1:length(D(simulationToPlot).receptionModes(sampleToPlot).SPLOS)
    if isempty(D(simulationToPlot).receptionModes(sampleToPlot).SPLOS)
        break
    end
    prn = D(simulationToPlot).receptionModes(sampleToPlot).SPLOS(s);
    colorIndex = find(prn==D(simulationToPlot).receivedReceptionModes(constellationRepetitionToPlot).list);
    lineColor = colorsToUse{colorIndex};
    currLOSvec = [results.satellitePositions_ENU(:, prn), results.satelliteArrivals_ENU(:, prn)]'; % maybe switch
    
    plot3(currLOSvec(:,1), currLOSvec(:,2), currLOSvec(:,3),...
        'color', lineColor, 'LineStyle', '-') 
end

for s = 1:1:results.noReflections %
    
    prn = results.reflections(s).prn;
    colorIndex =  find(prn==D(simulationToPlot).receivedReceptionModes(constellationRepetitionToPlot).list);
    lineColor = colorsToUse{colorIndex};
    currLOSvec = [results.satellitePositions_ENU(:, prn), results.satelliteArrivals_ENU(:, prn)]'; % maybe switch
    
    if results.satelliteLineOfSightFlags(prn) == 1
        plot3(currLOSvec(:,1), currLOSvec(:,2), currLOSvec(:,3),...
            'color', lineColor, 'LineStyle', '-') 
    end
    
    reflectionPoint_ENU = results.reflections(s).reflectionPoint_ENU;
    Bplot = [reflectionPoint_ENU, results.satellitePositions_ENU(:, prn)]';
    Rplot = [results.antennaPosition_ENU, reflectionPoint_ENU]';
    
    hold on
    scatter3(reflectionPoint_ENU(1), reflectionPoint_ENU(2), reflectionPoint_ENU(3),...
        'MarkerFaceColor', lineColor, 'Marker', '*')
    plot3(Bplot(:,1), Bplot(:,2), Bplot(:,3), 'color', lineColor, 'LineStyle', '-.')
    plot3(Rplot(:,1), Rplot(:,2), Rplot(:,3), 'color', lineColor, 'LineStyle', '-.')
    
end

clear prnLegend
h = zeros(1, length(D(simulationToPlot).receivedReceptionModes(constellationRepetitionToPlot).list));
for s = 1:1:length(D(simulationToPlot).receivedReceptionModes(constellationRepetitionToPlot).list)
    h(s) = plot3(NaN, NaN, NaN, 'color', colorsToUse{s});
    prnLegend{s} = ['PRN ' num2str(D(simulationToPlot).receivedReceptionModes(constellationRepetitionToPlot).list(s))];
end
legend(h, prnLegend);

disp('Signal paths plotted!')
%%
% Part 10: Simulate a receiver to compute pseudorange errors
results_folder = '../results'
samplingFactor = 0.1; % Sampling factor of the correlator output, ratio of the carrier frequency
correlatorSpacing = 0.5; % Spacing of the early late phase correlator, number of chips (between 0 and 1)
noChips_delay = 2; % Number of chips in each direction to find correlation peak
noChips_SNR = 30; % Number of chips in each direction to decide signal to noise ratio
phaseNoiseStandardDeviation = 0; % phase noise of receiver
amplitudeNoiseStandardDeviation = 0; % amplitude noise of receiver

constellation = 'GPS'; % The constellation for which to compute the pseudorange errors: GPS or GAL
frequencyBands = {'L1', 'L2'};
signalComponents = {'CA', 'CM'};


disp('Starting receiver simulation to compute pseudorange errors...')
noFrequencies = 1:1:length(frequencyBands);
for d = 1:1:length(D)

    clear measurements
    results = D(d).data.results;
    
    for c = 1:1:size(results,1)

        for s = 1:1:size(results,2) 

            clear buildingReflectionIndices
            for r = 1:1:length(results(c, s).reflections) % Exclude reflections on car
                buildingReflectionIndices(r) = ~strcmp(results(c, s).reflections(r).objectType, 'car');
            end
            reflections = results(c, s).reflections(buildingReflectionIndices);
            satellitesWithReflection = [reflections.prn];
            [C, IA, IC] = unique(satellitesWithReflection);
            
            modifiedObservations = zeros(0,6);
            for satellite = 1:1:length(IA) % Loop over all satellites with at least 1 reflection

                satellitePRN = C(satellite);
                newSatelliteIndices = IA(satellite):1:(IA(satellite) + sum(IC==satellite) - 1); 

                lineOfSightFlag = results(c, s).satelliteLineOfSightFlags(satellitePRN);
                
                noReflections = length(newSatelliteIndices);
                clear reflectionDelay powerLossRatio
                for r = 1:1:noReflections 
                    reflectionDelay(r, 1:2) = reflections(newSatelliteIndices(r)).reflectionDelay(1:noFrequencies);
                    powerLoss(r, 1:2) = reflections(newSatelliteIndices(r)).powerLoss(1:noFrequencies);
                end
                
                [pseudoRangeErrors, signalStrengthLossRatio] = computePseudorange(samplingFactor,...
                    correlatorSpacing, noChips_delay, noChips_SNR, ...
                    phaseNoiseStandardDeviation, amplitudeNoiseStandardDeviation, ...
                    constellation, frequencyBands, signalComponents, satellitePRN,...
                    lineOfSightFlag, noReflections, reflectionDelay, 10.^(-powerLoss/20), const);
                
                if lineOfSightFlag == 1
                    modeFlag = 2; % MP
                else
                    modeFlag = 3; % NLOS
                end
                
                satelliteElevation = results(c,s).satellitePositions_AER(2,satellitePRN);
                SNR_L1 = getSNRfromElevation(satelliteElevation);
                SNR_L2 = SNR_L1 - 3;
                
                modifiedSNR_L1 = 20*log10(10.^(SNR_L1/20) * signalStrengthLossRatio(1));
                modifiedSNR_L2 = 20*log10(10.^(SNR_L2/20) * signalStrengthLossRatio(2));
                           
                modifiedObservations(satellite,1:1:6) = [satellitePRN,...
                    pseudoRangeErrors(1), pseudoRangeErrors(2),...
                    modifiedSNR_L1, modifiedSNR_L2, modeFlag];
            end
                
            satelliteReceiverDistances = sqrt(sum(...
                [D(d).data.results(c, s).antennaPosition_ECEF(1) - D(d).data.satelliteTrajectories.positions_ECEFx(:, s),...
                D(d).data.results(c, s).antennaPosition_ECEF(2) - D(d).data.satelliteTrajectories.positions_ECEFy(:, s),...
                D(d).data.results(c, s).antennaPosition_ECEF(3) - D(d).data.satelliteTrajectories.positions_ECEFz(:, s)]...
                .^2,2)) ;
            
            
            % Generate measurements for open area (all satellites above 10
            % degree elevation are received, no pseudorange errors)
            satellitesAbove = find(D(d).data.results(c, s).satellitePositions_AER(2,:) >= 10 );
            
            satelliteElevations = results(c, s).satellitePositions_AER(2,:);
            SNRs_L1 = getSNRfromElevation(satelliteElevations);
            SNRs_L2 = SNRs_L1 - 3;
            
            measurements.open(c, s).pseudorange_L1 = ...
                satelliteReceiverDistances(satellitesAbove);
            measurements.open(c, s).pseudorange_L2 = ...
                satelliteReceiverDistances(satellitesAbove);
            
            measurements.open(c, s).satellitePositions = [...
                D(d).data.satelliteTrajectories.positions_ECEFx(satellitesAbove, s),...
                D(d).data.satelliteTrajectories.positions_ECEFy(satellitesAbove, s),...
                D(d).data.satelliteTrajectories.positions_ECEFz(satellitesAbove, s)];
            
            measurements.open(c, s).PRNs = satellitesAbove;
            measurements.open(c, s).receptionModes = ones(1,length(satellitesAbove))';
            % receptionModes 1 = SPLOS
            % receptionModes 2 = MP
            % receptionModes 3 = NLOS
            measurements.open(c, s).SNRs_L1 = SNRs_L1(satellitesAbove)';
            measurements.open(c, s).SNRs_L2 = SNRs_L2(satellitesAbove)';
            measurements.open(c, s).elevations = satelliteElevations(satellitesAbove)';

            
            % Generate measurements for the urban area
            PRNs_SPLOS = intersect(find(D(d).data.results(c, s).satelliteLineOfSightFlags), ...
                setdiff(1:1:31,satellitesWithReflection)); % 31 or 55?

            pseudoranges_L1_SPLOS = satelliteReceiverDistances(PRNs_SPLOS);
            pseudoranges_L2_SPLOS = satelliteReceiverDistances(PRNs_SPLOS);
            SNRs_L1_SPLOS = SNRs_L1(PRNs_SPLOS)';
            SNRs_L2_SPLOS = SNRs_L2(PRNs_SPLOS)';
            elevetations_SPLOS = satelliteElevations(PRNs_SPLOS)';
            satellitePositions_SPLOS = [D(d).data.satelliteTrajectories.positions_ECEFx(PRNs_SPLOS, s),...
                D(d).data.satelliteTrajectories.positions_ECEFy(PRNs_SPLOS, s),...
                D(d).data.satelliteTrajectories.positions_ECEFz(PRNs_SPLOS, s)];
            receptionModes_SPLOS = ones(1, length(pseudoranges_L1_SPLOS));
            
            PRNs_reflected = modifiedObservations(:,1)';
            pseudoranges_L1_reflected = satelliteReceiverDistances(PRNs_reflected) + modifiedObservations(:,2);
            pseudoranges_L2_reflected = satelliteReceiverDistances(PRNs_reflected) + modifiedObservations(:,3);
            SNR_L1_reflected = modifiedObservations(:,4);
            SNR_L2_reflected = modifiedObservations(:,5);
            satellitePositions_reflected = [D(d).data.satelliteTrajectories.positions_ECEFx(PRNs_reflected, s),...
                D(d).data.satelliteTrajectories.positions_ECEFy(PRNs_reflected, s),...
                D(d).data.satelliteTrajectories.positions_ECEFz(PRNs_reflected, s)];
            receptionModes_reflected = modifiedObservations(:,6)';
            elevations_reflected = satelliteElevations(PRNs_reflected)';
            
            PRNs_urban = [PRNs_SPLOS, PRNs_reflected];
            pseudoranges_L1_urban = [pseudoranges_L1_SPLOS; pseudoranges_L1_reflected];
            pseudoranges_L2_urban = [pseudoranges_L2_SPLOS; pseudoranges_L2_reflected];
            SNRs_L1_urban = [SNRs_L1_SPLOS; SNR_L1_reflected];
            SNRs_L2_urban = [SNRs_L2_SPLOS; SNR_L2_reflected];
            satellitePositions_urban = [satellitePositions_SPLOS; satellitePositions_reflected];
            receptionModes_urban = [receptionModes_SPLOS, receptionModes_reflected];
            elevations_urban = [elevetations_SPLOS; elevations_reflected];

            [measurements.urban(c, s).PRNs, order] = sort(PRNs_urban);
            measurements.urban(c, s).pseudorange_L1 = pseudoranges_L1_urban(order);
            measurements.urban(c, s).pseudorange_L2 = pseudoranges_L2_urban(order);
            measurements.urban(c, s).satellitePositions = satellitePositions_urban(order,:);
            measurements.urban(c, s).receptionModes = receptionModes_urban(order)';
            
            measurements.urban(c, s).SNRs_L1 = SNRs_L1_urban(order);
            measurements.urban(c, s).SNRs_L2 = SNRs_L2_urban(order);
            measurements.urban(c, s).elevations = elevations_urban(order);
            
        end
    end
    
    results_fileName = ['receiver' num2str(d)];
    save([results_folder '/' results_fileName '.mat'], 'measurements')
end
disp('New reflection-influenced pseudoranges computed and saved!')
%% 
% Part 11: Compute and plot pseudorange error distribution
init

tempData_folder = '../Data_temp/ucSim/rec/'
results_folder = '../results'

filesNumberToAnalyze = [1]

font = 'latex' % Interpreter for the text in the plots
fontSize = 30 % Font size of the text in the plots
plotLineWidth = 2

set(0, 'DefaultTextInterpreter', font)
set(0, 'DefaultLegendInterpreter', font)
set(0, 'DefaultAxesTickLabelInterpreter', font)
set(0, 'defaultTextFontSize', fontSize)
set(0, 'defaultAxesFontSize', fontSize)
set(0, 'defaultlinelinewidth', plotLineWidth);

pseudorangeErrors_frequency1 = [];
pseudorangeErrors_frequency2 = [];
elevations = [];
receptionModes  = [];
for d = filesNumberToAnalyze 
    
    tempData_fileName = ['receiver' num2str(d)];
    load([results_folder '/' tempData_fileName '.mat'])

    for c = 1:1:size(measurements.open, 1)    
        for s = 1:1:size(measurements.open, 2)
            openResults = measurements.open(c, s);
            urbanResults = measurements.urban(c, s);
            tooLowSatellites = find(urbanResults.elevations < 10);
            if ~isempty(tooLowSatellites)
                urbanResults.pseudorange_L1(tooLowSatellites) = [];
                urbanResults.pseudorange_L2(tooLowSatellites) = [];
                urbanResults.elevations(tooLowSatellites) = [];
                urbanResults.recMode(tooLowSatellites) = [];
            end
            
            if ~isempty(setdiff(openResults.PRNs, urbanResults.PRNs))
                [~, toRemove] = ismember(setdiff(openResults.PRNs, urbanResults.PRNs), openResults.PRNs );
                openResults.pseudorange_L1(toRemove) = [];
                openResults.pseudorange_L2(toRemove) = [];
                openResults.elevations(toRemove) = [];
                openResults.receptionModes(toRemove) = [];
            end
            pseudorangeErrors_frequency1 = [pseudorangeErrors_frequency1; urbanResults.pseudorange_L1 - openResults.pseudorange_L1];
            pseudorangeErrors_frequency2 = [pseudorangeErrors_frequency2; urbanResults.pseudorange_L2 - openResults.pseudorange_L2];
            elevations = [elevations; openResults.elevations];
            receptionModes = [receptionModes; urbanResults.receptionModes];
        end
    end
end

close all
figure('units', 'normalized', 'outerposition', [0 0 1 1])
histogram(pseudorangeErrors_frequency1, 'normalization', 'probability');
ylabel('Probability')
xlabel('Pseudorange error on frequency 1 (m)')

figure('units', 'normalized', 'outerposition', [0 0 1 1])
histogram(pseudorangeErrors_frequency2, 'normalization', 'probability');
ylabel('Probability')
xlabel('Pseudorange error on frequency 2 (m)')