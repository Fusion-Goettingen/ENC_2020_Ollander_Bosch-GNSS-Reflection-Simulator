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
function crossingPoint = planeLineIntersect(normal, planePoint, planeLimits,...
    linePoint1, linePoint2)
% PLANELINEINTERSECT finds the intersection point between a line and a plane (if possible)
% Input:
%        normal,            3x1, normal vector of the plane
%        planePoint,        3x1, a point on the plane
%        planeLimits,       3x2, limits of the plane
%        linePoint1,        3x1, first point on the line
%        linePoint2,        3x1, second point on the line
%
% Output:
%        crossingPoint,     3x1, the identified crossing point (or NaN if there is none)
%
% Written by Simon Ollander

% Intersect plane and line
crossingPoint = linePoint2 + sum((planePoint-linePoint2).* normal) /...
    sum((linePoint1 - linePoint2).* normal) * (linePoint1 - linePoint2);

sensitivity = 1e-3; % 1 mm

% Check that crosssPoint is within plane limits
crossPointWithinPlane = NaN*ones(3,1);
for d = 1:1:3
    crossPointWithinPlane(d) = crossingPoint(d) >= planeLimits(d,1) - sensitivity &&...
        crossingPoint(d) <= planeLimits(d,2) + sensitivity;
    if ~crossPointWithinPlane(d)
        crossingPoint = NaN*ones(3,1);
        return
    end
end

% Check that cross point lies between satellite and antenna
v = (linePoint2 - linePoint1);
t = (crossingPoint - linePoint1)' / v';
crossIsPointBetweenAntennaAndSatellite = (t >= 0 && t <= 1);

if ~crossIsPointBetweenAntennaAndSatellite
    crossingPoint = NaN*ones(3,1);
end

end