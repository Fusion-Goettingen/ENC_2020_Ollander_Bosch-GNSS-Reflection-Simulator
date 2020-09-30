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
function outFrequency = computeDoppler(receiverPosition, receiverVelocity, ...
    emitterPosition, emitterVelocity, inFrequency, const)
% COMPUTEDOPPLER Computes the output frequency after a Doppler shift, given
% 3D receiver and emitter positions and velocities
% Input:
%        receiverPosition,	3x1, position of the receiver (meters)
%        receiverVelocity,	3x1, velocity of the receiver (meters/second)
%        emitterPosition,	3x1, position of the emitter (meters)
%        emitterVelocity,	3x1, velocity of the emitter (meters/second)
%        inFrequency,       1x1, carrier frequency of the incoming signal (Hz)
%        const,             struct containing natural constants (speed of light)
%
% Output:
%        outFrequency,      1x1, carrier frequency of the Doppler-shifted signal (Hz)       
%
% Written by Simon Ollander
%
% For more details and visualization:
% https://isaacphysics.org/concepts/cp_doppler_effect
% https://en.wikipedia.org/wiki/File:DopplerSatScheme.png

antennaSatelliteVector = emitterPosition - receiverPosition;

outFrequency = inFrequency .* (const.c * norm(antennaSatelliteVector) - dot(receiverVelocity,antennaSatelliteVector)) /...
    (const.c * norm(antennaSatelliteVector)  - dot(emitterVelocity,antennaSatelliteVector));

end