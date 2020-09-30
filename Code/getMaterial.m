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
function material = getMaterial(reflectorType)
% GETMATERIAL returns the material of a given reflector type
% Input:
%        reflectorType, string defining the type of the reflector: 'build', 'ground', 'car'
%
% Output:
%        material,      string defining the material of the reflector: 'concrete', 'dryGround', 'steel', 'glass'
%
% Written by Simon Ollander
if strcmp(reflectorType, 'build')
    %Glass 'glass' could also be incorporated here, or we could randomly select concrete or glass
    material = 'concrete';
elseif strcmp(reflectorType, 'ground')
    material = 'dryGround';
elseif strcmp(reflectorType, 'car')
    material = 'steel';
end

end