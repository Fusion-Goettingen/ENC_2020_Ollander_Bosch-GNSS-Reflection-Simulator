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
function reflectionCoefficient = computeReflectionCoefficient(const, material, carrierFrequencies,...
    reflectionAngle, rejectionLHCP, elevation) 
% COMPUTEREFLECTIONCOEFFICIENT outputs the reflection coefficient for one reflection
% Input:
%        const,                     struct containing natural constants
%        material,                  string defining the material of the reflection surface ('concrete')
%        carrierFrequencies,        1xN, vector defining the carrier frequencies (Hz)
%        reflectionAngle,           1x1, reflection angle (radians)
%        rejectionLHCP,             1x1, left hand circual polarization rejection of the antenna (dB)
%        elevation,                 1x1, elevation angle of the satellite (degrees)
%
% Output:
%        reflectionCoefficient,    	struct containing reflection coeffients (horizontal, vertical, copolar, crosspolar, effective)
%
% Written by Simon Ollander

relativePermittivity = const.relativePermittivity.(material);
dielectricConstant = const.conductivity.(material);
wavelengths = const.c ./ carrierFrequencies; 

epsil = relativePermittivity - 1j*60*wavelengths*dielectricConstant; % Complex dielectric constants
reflectionCoefficient.Horizontal = (sin(reflectionAngle) - sqrt(epsil - cos(reflectionAngle).^2)) ./ (sin(reflectionAngle) + sqrt(epsil - cos(reflectionAngle).^2)); % Horizontal reflection coefficient
reflectionCoefficient.Vertical = (epsil*sin(reflectionAngle) - sqrt(epsil - cos(reflectionAngle).^2)) ./ (epsil*sin(reflectionAngle) + sqrt(epsil - cos(reflectionAngle).^2)); % Vertical reflection coefficient
reflectionCoefficient.Copolar = (reflectionCoefficient.Horizontal + reflectionCoefficient.Vertical) / 2;
reflectionCoefficient.Crosspolar = (reflectionCoefficient.Horizontal - reflectionCoefficient.Vertical) / 2 ;
reflectionCoefficient.Effective = (abs(reflectionCoefficient.Copolar) + 10^(-rejectionLHCP*sind(elevation)/20)*abs(reflectionCoefficient.Crosspolar))*exp(-j*pi); % Effective reflection coefficient

end