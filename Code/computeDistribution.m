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
function D = computeDistribution(D)
% COMPUTEDISTRIBUTION goes through all reflections and summarizes the
% results according to the reception mode
% Input:
%        D,    struct containing all reflections
%
% Output:
%        D,    struct containing all reflections, including reception modes and distributions
%
% Written by Simon Ollander
reflectionParameters = {'reflectionDelay', 'powerLoss','dopplerDifference', 'phaseDifference'};

for d = 1:1:length(D)
    clear reflectionDistribution
    
    addedReflections = 0;
    for c = 1:1:size(D(d).data.results,1)
        for s = 1:1:size(D(d).data.results,2)
            if ~isempty(D(d).data.results(c, s).reflections)
                results = D(d).data.results(c, s);
                for r = 1:1:length(D(d).data.results(c, s).reflections)
                    reflection = results.reflections(r);
                    addedReflections = addedReflections + 1;
                    
                    reflectionDistribution.LOS(addedReflections) =...
                        results.satelliteLineOfSightFlags(reflection.prn);
                    reflectionDistribution.planeType{addedReflections} =...
                        reflection.planeType;
                    reflectionDistribution.objectType{addedReflections} =...
                        reflection.objectType;
                    reflectionDistribution.azimuth(addedReflections) = ...
                        results.satellitePositions_AER(1, reflection.prn);
                    reflectionDistribution.elevevation(addedReflections) = ...
                        results.satellitePositions_AER(2, reflection.prn);
                    
                    for field = 1:1:length(reflectionParameters)
                        newData = reflection.(reflectionParameters{field});
                        reflectionDistribution.(reflectionParameters{field})(1:1:3, addedReflections) = newData';
                    end
                end
            end
        end
    end
    D(d).distrubution.reflections = reflectionDistribution;
end

for d = 1:1:length(D)
    results = D(d).data.results;
    clear receptionModes receivedPRNs
    
    % List of PRNs that are received during the entire simulation (for color assignment in plot)
    for c = 1:1:size(D(d).data.results, 1)
        receivedPRNs(c).list = zeros(1,0) ;
    end
    
    D(d).numberOfReflections = 0;
    for c = 1:1:size(results, 1)
        for s = 1:1:size(results, 2)
            noReflectionsToAdd = 0;
            if ~isempty(results(c, s).reflections)
                noReflectionsToAdd = length(find(strcmp({results(c, s).reflections.objectType}, 'build')...
                    & strcmp({results(c, s).reflections.planeType} , 'wall')));
            end
            D(d).numberOfReflections = D(d).numberOfReflections  + noReflectionsToAdd;
            
            satellitesBelow = find(results(c, s).satellitePositions_AER(2,:) < 0);
            satellitesLOS = find(results(c, s).satelliteLineOfSightFlags);
            satellitesNLOS = find(~results(c, s).satelliteLineOfSightFlags);
            
            receivedPRNs(c).list =...
                union(receivedPRNs(c).list, satellitesLOS);
            
            reflectedPRNs = [];
            for r = 1:1:results(c, s).noReflections
                if strcmp(results(c,s).reflections(r).objectType, 'build') ...
                        && strcmp(results(c,s).reflections(r).planeType, 'wall')
                    reflectedPRNs = [reflectedPRNs, results(c,s).reflections(r).prn];
                end
            end
            
            receivedPRNs(c).list = union(receivedPRNs(c).list, reflectedPRNs);
            
            receptionModes(c, s).SPLOS = setdiff(satellitesLOS, reflectedPRNs);
            receptionModes(c, s).MP = intersect(reflectedPRNs, satellitesLOS);
            receptionModes(c, s).NLOS = intersect(reflectedPRNs, satellitesNLOS);
            receptionModes(c, s).Blocked = setdiff(setdiff(satellitesNLOS, satellitesBelow), reflectedPRNs);
            receptionModes(c, s).tooLow = satellitesBelow;
            
        end
    end
    
    D(d).receivedReceptionModes = receivedPRNs;
    D(d).receptionModes = receptionModes;
end

% Compute distribution of reception modes
receptionModeList = {'SPLOS', 'MP', 'NLOS', 'Blocked', 'tooLow'};
for d = 1:1:length(D)
    receptionModeDistribution = zeros(1, length(receptionModeList));
    for c = 1:1:size(D(d).data.results, 1)
        for s = 1:1:size(D(d).data.results, 2)
            for p = 1:1:length(receptionModeList)
                receptionModeDistribution(p) = receptionModeDistribution(p) + ...
                    length(D(d).receptionModes(c, s).(receptionModeList{p}));
            end
        end
    end
    D(d).receptionModeDistribution = receptionModeDistribution /...
        (size(D(d).data.results, 1)*size(D(d).data.results, 2));
    
end

end