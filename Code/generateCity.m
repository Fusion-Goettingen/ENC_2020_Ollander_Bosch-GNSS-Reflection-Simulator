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
function [planes, blockCorners] = generateCity(geometry)
% GENERATECITY generates an urban canyon environment, consisting of planes that make up buildings
% Input:
%        geometry,    struct defining the geometrical configuration of the urban canyon
%
% Output:
%        planes,        list of the planes in the urban canyon
%        blockCorners,	struct defining the corners of each city block
%
% Written by Simon Ollander
lastCorner = 0;
blockCorners.x = NaN*ones(3, 1);
for b = 1:1:geometry.block.noBlocks.x
    blockCorners.x(b) = lastCorner + geometry.road.width + geometry.block.size.x;
    lastCorner = blockCorners.x(b);
end
blockCorners.x = blockCorners.x - geometry.block.size.x;

lastCorner = 0;
blockCorners.y = NaN*ones(3,1);
for b = 1:1:geometry.block.noBlocks.y
    blockCorners.y(b) = lastCorner + geometry.road.width + geometry.block.size.y;
    lastCorner = blockCorners.y(b);
end
blockCorners.y = blockCorners.y - geometry.block.size.y;

geometry.block.sizes = [geometry.block.size.x, geometry.block.size.y, 0]';

planes = zeros(0, 1);
for bx = 1:1:geometry.block.noBlocks.x
    xFirst = bx == 1;
    xLast = bx == geometry.block.noBlocks.x;
    
    for by = 1:1:geometry.block.noBlocks.y
        
        if (geometry.skipCentralBlock == 1) && (bx == ceil(geometry.block.noBlocks.x/2)) && (by == ceil(geometry.block.noBlocks.y/2))
            % Skip middle block to create urban square / town square environment
        else       
            yFirst = by == 1;
            yLast = by == geometry.block.noBlocks.y;
            
            block.startPoint = geometry.centerShift + [blockCorners.x(bx), blockCorners.y(by), 0]';
            sideSkips(1).toSkip = [];
            sideSkips(2).toSkip = [];
            
            if yFirst
                sideSkips(1).toSkip = 1;
            elseif yLast
                sideSkips(1).toSkip = 2;
            end
            
            if xFirst
                sideSkips(2).toSkip = 1;
            elseif xLast
                sideSkips(2).toSkip = 2;
            end
            
            for currDimension = 1:2
                otherDim = setdiff(1:2, currDimension);
                currSkip =  sideSkips(currDimension).toSkip;
                for currSide = setdiff(1:2,currSkip)
                    
                    if currSide == 1
                        build.startPoint = [0, 0, 0]';
                    elseif currSide == 2 && currDimension == 1
                        build.startPoint = [0, geometry.block.sizes(2)-1*geometry.build.depth, 0]';
                    elseif currSide == 2 && currDimension == 2
                        build.startPoint = [geometry.block.sizes(1)-1*geometry.build.depth, 0, 0]';
                    end
                    
                    build.startPoint = build.startPoint + block.startPoint;
                    
                    while build.startPoint(currDimension) < geometry.block.sizes(currDimension) +  block.startPoint(currDimension)
                        if geometry.build.deviationWidth == 0
                            newWidth = geometry.build.meanWidth;
                        else
                            newWidth =  random('Rician', geometry.build.meanWidth, geometry.build.deviationWidth);
                        end
                        
                        if newWidth > geometry.block.sizes(currDimension) - build.startPoint(currDimension) +  block.startPoint(currDimension)
                            newWidth = geometry.block.sizes(currDimension) - build.startPoint(currDimension) +  block.startPoint(currDimension);
                        end
                        build.size = NaN*ones(3,1);
                        build.size(currDimension) = newWidth;
                        
                        if geometry.build.deviationHeight == 0
                            newHeight = geometry.build.meanHeight;
                        else
                            newHeight = random('Rician', geometry.build.meanHeight, geometry.build.deviationHeight);
                        end
                        build.size(3) = newHeight;
                        
                        build.size(otherDim) = geometry.build.depth;
                        
                        build.upper = build.startPoint + build.size;
                        
                        depthOffset_pos = 1*geometry.build.maxDepthOffset*rand; % Uniform distribution
                        buildStartPointTemp = build.startPoint;
                        if currDimension == 1 && currSide == 1
                            buildStartPointTemp(2) = buildStartPointTemp(2) + depthOffset_pos;
                        elseif currDimension == 1 && currSide == 2
                            build.upper(2) = build.upper(2) - depthOffset_pos;
                        elseif currDimension == 2 && currSide == 1
                            buildStartPointTemp(1) = buildStartPointTemp(1) + depthOffset_pos;
                        elseif currDimension == 2 && currSide == 2
                            build.upper(1) = build.upper(1) - depthOffset_pos;
                        end
                        
                        newPlanes = generateBlock(buildStartPointTemp, build.upper, 'build');
                        planes = [planes, newPlanes];
                        
                        newLimits = [0, 0, 0]';
                        newLimits(currDimension) = build.size(currDimension);
                        build.startPoint = build.startPoint + newLimits;
                    end
                end
            end
        end
    end
end

upperPlotLimit = max([blockCorners.x(end) + geometry.block.size.x + geometry.road.width...
    blockCorners.y(end) + geometry.block.size.y + geometry.road.width]);

% Add the ground as reflection plane
extraPlaneLimits = 1e4;
groundPlane.normal = [0, 0, 1]';
groundPlane.limits = [-extraPlaneLimits, upperPlotLimit+extraPlaneLimits;...
    -extraPlaneLimits, upperPlotLimit+extraPlaneLimits;...
    -extraPlaneLimits, extraPlaneLimits];
groundPlane.points = [-extraPlaneLimits, -extraPlaneLimits, 0 ;...
    upperPlotLimit+extraPlaneLimits, -extraPlaneLimits, 0;...
    upperPlotLimit+extraPlaneLimits, upperPlotLimit+extraPlaneLimits, 0;...
    -extraPlaneLimits, upperPlotLimit+extraPlaneLimits, 0]';
groundPlane.planeType = 'floor';
groundPlane.objectType = 'ground';
planes(end+1) = groundPlane;

end