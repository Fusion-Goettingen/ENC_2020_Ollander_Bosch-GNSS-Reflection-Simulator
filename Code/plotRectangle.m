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
function plotRectangle(point1, point2, point3, point4, color)
% PLOTRECTANGLE plots a rectangle in 3D
% Input:
%        point1,    3x1, first corner of the rectangle
%        point2,    3x1, second corner of the rectangle
%        point3,    3x1, third corner of the rectangle
%        point4,    3x1, fourth corner of the rectangle
%        color,     char or 3x1 vector defining the color of the rectangle
%
% Written by Simon Ollander

x = [point1(1), point2(1), point3(1), point4(1)];
y = [point1(2), point2(2), point3(2), point4(2)];
z = [point1(3), point2(3), point3(3), point4(3)];
fill3(x, y, z, color);
hold on

end