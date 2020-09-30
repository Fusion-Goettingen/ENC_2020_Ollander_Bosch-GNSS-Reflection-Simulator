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
function SNR = getSNRfromElevation(elevation)
% GETSNRFROMELEVATION returns the signal-to-noise ratio based on the satellite elevation
% Input:
%        elevation,    1x1, the elevation of the satellite (degrees)
%
% Output:
%        SNR,           1x1, the elevation-based signal-to-noise ratio (dB)
%
% Written by Simon Ollander
%
% Source of the model: https://www.sciencedirect.com/science/article/pii/S1110982316300412
SNR = 3.199e-5 * elevation.^3 -0.0081*elevation.^2 +0.6613*elevation + 31.38;

end