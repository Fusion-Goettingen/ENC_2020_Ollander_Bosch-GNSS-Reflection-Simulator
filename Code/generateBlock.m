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
function newPlanes = generateBlock(underCorner, upperCorner, blockType)
% GENERATEBLOCK uses two corners to generate a 3D block consisting of 4 walls and one roof
% Input:
%        underCorner,    3x1, lower corner of the block to be generated (meters)
%        upperCorner,    3x1, upper corner of the block to be generated (meters)
%        blockType,      string defining the type of block ('car', 'building'...)
%
% Output:
%        newPlanes,      struct contaning a list with the new planes representing the block (4 walls and 1 roof)
%
% Written by Simon Ollander
%
% For visualization, see http://www.matrixlab-examples.com/3d-polygon.html
p1 = [underCorner(1), underCorner(2), underCorner(3)]';
p2 = [upperCorner(1), underCorner(2), underCorner(3)]';
p3 = [upperCorner(1), upperCorner(2), underCorner(3)]';
p4 = [underCorner(1), upperCorner(2), underCorner(3)]';
p5 = [underCorner(1), underCorner(2), upperCorner(3)]';
p6 = [upperCorner(1), underCorner(2), upperCorner(3)]';
p7 = [upperCorner(1), upperCorner(2), upperCorner(3)]';
p8 = [underCorner(1), upperCorner(2), upperCorner(3)]';

myWalls(1).points = [p5, p6, p7, p8]; % roof
%walls 1-4:
myWalls(2).points = [p2, p6, p5, p1]; % opposite myWalls(4)
myWalls(3).points = [p2, p6, p7, p3]; % opposite myWalls(5)
myWalls(4).points = [p3, p4, p8, p7];
myWalls(5).points = [p1, p4, p8, p5];

myWalls(6).points = [p1, p2, p3, p4]; % floor

oppositeTable(1) = 6;
oppositeTable(2) = 4;
oppositeTable(3) = 5;
oppositeTable(4) = 2;
oppositeTable(5) = 3;
oppositeTable(6) = 1;

w = 0;
for m = 1:1:length(myWalls)
    
    wp1 = myWalls(m).points(:, 1);
    wp2 = myWalls(m).points(:, 2);
    wp3 = myWalls(m).points(:, 3);
    wp4 = myWalls(m).points(:, 4);
    normal = cross(wp1-wp2, wp1-wp3);
    normal = normal/norm(normal);
    
    % Force normal outwards:
    oppositePoint = myWalls(oppositeTable(m)).points(:,1);
    
    % Vector from opposite wall to current wall, should align with normal:
    oppositeVector = wp1 - oppositePoint;
    
    if abs(acos( sum(oppositeVector.* normal) / (norm(oppositeVector) * norm(normal)))) > pi/2
        normal = -1*normal; % invert direction of normal
    end
    
    allwp = [wp1, wp2, wp3, wp4];
    limits = [min(allwp(1,:)), max(allwp(1,:));...
        min(allwp(2,:)), max(allwp(2,:));...
        min(allwp(3,:)), max(allwp(3,:))];
    
    % Add new planes to list: (ignore floors)
    if m < 6
        w = w + 1;
        newPlanes(w).normal = normal;
        newPlanes(w).limits = limits;
        newPlanes(w).points = [wp1, wp2, wp3, wp4];
        
        % Specifiy wall type
        if m == 1
            newPlanes(w).planeType = 'roof';
        elseif m == 6
            newPlanes(w).planeType = 'floor';
        else
            newPlanes(w).planeType = 'wall';
        end
        newPlanes(w).objectType = blockType;
    end
end


end