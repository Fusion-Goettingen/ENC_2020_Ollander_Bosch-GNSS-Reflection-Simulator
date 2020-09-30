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
function [reflectionPoint, mirrorPoint, projectionPoint] =....
    computeReflection(normal, planePoint, planeLimits, antennaPosition,...
    satellitePosition)
% COMPUTEREFLECTION Uses orthogonal projection to compute if a reflection between two points is possible in a plane
% Input:
%        normal,            3x1, normal vector of the reflection plane
%        planePoint,        3x1, a point of the reflection plane
%        planeLimits,       3x2, under and upper limits of the plane
%        antennaPosition,   3x1, position of the antenna
%        satellitePosition, 3x1, position of the satellite
%
% Output:
%        reflectionPoint,   3x1, point of the reflection
%        mirrorPoint,       3x1, mirror point of the receiver in the plane
%        projectionPoint,   3x1, projection point of the satellite in the plane
%
% Written by Simon Ollander

projectionPoint = antennaPosition + ...
    sum((planePoint-antennaPosition ).* normal) ...
    / sum(normal .* normal) * normal;

mirrorPoint = antennaPosition + 2*(projectionPoint-antennaPosition);

reflectionPoint = planeLineIntersect(normal, planePoint, planeLimits, ...
    mirrorPoint, satellitePosition);

end