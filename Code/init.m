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
% This script initializes the program and loads all natural constants into the variable const
clear all; clc; close all;

addpath ../RangingCodes

% Define constants
const.c = 299792458; % Speed of light, m/s
const.eps0 = 8.854187817e-12; % Permeability in vacuum, F/m
const.R_e = 6371008.7714; % Earth radius, m

% GPS orbital information
const.orbits.GPS.semiMajorAxis = [5153.653320^2, 5153.606445^2, 5153.542969^2, 5153.562012^2, 5153.740723^2, 5153.577637^2, 5153.717285^2, 5153.704102^2, 5153.624512^2, 5153.617188^2, 5153.59912^2, 5153.634766^2, 5153.541016^2, 5153.624512^2, 5153.543457^2, 5153.713867^2, 5153.660645^2, 5153.585938^2, 5153.533691^2, 5153.751465^2, 5153.631348^2, 5153.538086^2, 5153.669922^2, 5153.639160^2, 5153.549316^2, 5153.604980^2, 5153.591309^2, 5153.683105^2, 5153.708008^2, 5153.531250^2, 5153.600586^2];
const.orbits.GPS.meanAnomaly = [-0.1208933778e+001, -0.7379312701e+000, -0.2476050596e+001 0.1729284331e+001 -0.5955965034e+000 -0.1863089325e+001 0.1717615066e+001 0.1041308928e+001 0.3003120927e+001 -0.1772399827e+001 -0.2456528294e+001 -0.8542115828e+000 0.4893275126e+000 0.1755864967e+000 0.1876034910e+001 0.8991835103e+000 0.2479146271e+001 -0.2239168159e+001 0.5016783799e-001 -0.2845982272e+001 0.4089099929e+000 -0.3331816757e+000 -0.8946785652e+000 -0.2992513015e+001  0.2928989882e+001 0.9610869005e+000 0.1794198757e+001 0.2993612193e+001 -0.1855898041e+001 -0.2008106311e+001 0.1461498698e+001];
const.orbits.GPS.perigeeArgument = [0.563527841, -2.074941210 0.718299986 0.536422146 -1.700786593 -2.620313325 -1.011684673 1.869199782 -2.911753571  1.589303341 0.783939084  1.923932861 -1.944528126  0.517680691 0.382938305 -1.917666234 -1.854486524 0.879175098 1.466748537 -1.748448979 -1.854394021 -2.525418487 0.394514692 0.771921528 -0.184188174 0.328577112 -1.607607372 0.029954945 -3.122813373 -0.369073306 -2.439273631];
const.orbits.GPS.orbitalInclination = [0.9651879023, 0.9440357453 0.9595553166 0.9465644168 0.9649961547 0.9644988093 0.9649302414 0.9547136897 0.9597590484 0.8981062034  0.9901450506 0.9709223539 0.9634382054 0.9304156737 0.9906543802 0.9780709440 0.9252444807 0.9752067143 0.9265148085 0.9373066030 0.9232011703 0.9456176630 0.9475830759 0.9782207469 0.9601005988 0.9730315775 0.9893960366 0.9790955953 0.9502675423 0.9692924993 0.9592676952];
const.orbits.GPS.rightAscensionAscendingNode = [0.1489908052e+001, 0.1441165588e+001 0.2532003819e+001 0.2518117098e+001 0.1481505987e+001 -0.1634831037e+001 0.4283431626e+000 -0.2713451725e+001 0.2528123551e+001 0.1102367430e+001  -0.5551591071e+000 -0.2590521291e+001 -0.2628516527e+001 -0.2769660736e+001 -0.5363502406e+000 0.4829410443e+000 0.2493819082e+001 0.5313011368e+000 0.2442215384e+001 0.1452076850e+001 0.2494154640e+001 -0.2709604040e+001 -0.1684579427e+001 -0.6078234089e+000 -0.6194110309e+000 0.4322024576e+000 -0.5315595467e+000  0.4928152966e+000 -0.1592779889e+001 -0.1624817468e+001 -0.2711327521e+001];
const.orbits.GPS.eccentricity = [0.5968570709e-002, 0.1603031158e-001 0.2923011780e-003 0.4848480225e-002 0.2102851868e-003 0.1005315781e-001 0.1864910126e-002 0.6208419800e-003 0.1812458038e-002  0.1681089401e-001 0.6042003632e-002  0.4010200500e-002 0.8724212646e-002 0.8650302887e-002 0.8711338043e-002 0.1133203506e-001 0.1726150513e-001 0.1039505005e-001 0.4734039307e-002 0.2335453033e-001 0.7355213165E-002 0.1113939285e-001 0.4827976227e-002 0.5622386932e-002 0.1107215881E-002 0.3872394562e-002 0.2001571655e-001 0.5064010620e-003 0.2400398254e-002 0.8392333984e-002 0.6222724915e-003];
const.orbits.GPS.prnOrder = 1:1:length(const.orbits.GPS.semiMajorAxis);

% Galileo orbital information: https://www.gsc-europa.eu/system-service-status/orbital-and-technical-parameters
const.orbits.GAL.semiMajorAxis = 29599.8*1e3*ones(1,24);
const.orbits.GAL.meanAnomaly = pi/180*[15.153, 60.153, 345.153, 30.153, 150.153, 285.153, 135.153, 0.153, 120.153, 255.153, 225.153, 45.153, 75.153, 165.153, 300.153, 210.153 , 70.153, 90.153, 315.153, 180.153, 330.153, 195.153, 240.153, 105.153];
const.orbits.GAL.perigeeArgument = zeros(1,24);
const.orbits.GAL.orbitalInclination = 56.0*pi/180*ones(1,24);
const.orbits.GAL.rightAscensionAscendingNode = pi/180*[77.632, 77.632, 197.632, 197.632, 77.632, 77.632, 317.632, 317.632, 197.632, 197.632, 317.632, 317.632, 197.632, 197.632, 197.632, 197.632, 317.632, 317.632, 317.632, 317.632, 77.632, 77.632, 77.632, 77.632];
const.orbits.GAL.eccentricity = zeros(1,24);
const.orbits.GAL.prnOrder = [11, 12, 19, 20, 26, 22, 24, 30, 8, 9, 1, 2, 7, 3, 4, 5, 21, 25, 27, 31, 36, 13, 15, 33];

const.GPS.L1.CA.f_carrier = 1.57542e9; % Hz
const.GPS.L1.CA.T_code = 1e-3; % s
const.GPS.L1.CA.L_code = 1023; % bits
const.GPS.L1.CA.f_chip = 1.023e6; % bits/s
const.GPS.L1.CA.T_data = 50e-3; % s

const.GPS.L2.CL.f_carrier = 1.2276e9;
const.GPS.L2.CL.T_code = 1.5;
const.GPS.L2.CL.L_code = 767250;
const.GPS.L2.CL.f_chip = 0.5115e6;

const.GPS.L2.CM.f_carrier = 1.2276e9;
const.GPS.L2.CM.T_code = 20e-3;
const.GPS.L2.CM.L_code = 10230;
const.GPS.L2.CM.f_chip = 0.5115e6;

const.GPS.L5.SOL.f_carrier = 1176.45e6; % Secondary codes?????
const.GPS.L5.SOL.T_code = 1e-3;
const.GPS.L5.SOL.L_code = 10230;
const.GPS.L5.SOL.f_chip = 10.23e6;

const.GAL.E1.OS.f_carrier = 1.57542e9;
const.GAL.E1.OS.T_code = 4e-3;
const.GAL.E1.OS.L_code = 4092;
const.GAL.E1.OS.f_chip = 1.023e6;
const.GAL.E1.OS.f_subcarrier = 1.023e6;
const.GAL.E1.OS.m = 6;
const.GAL.E1.OS.n = 1;
const.GAL.E1.OS.alpha = sqrt(10/11);
const.GAL.E1.OS.beta = sqrt(1/11);

const.GAL.E5.OS.f_carrier = 1191.795e6;
const.GAL.E5.OS.T_code = 1e-3;
const.GAL.E5.OS.L_code = 10230;
const.GAL.E5.OS.f_chip = 15.345e6;
const.GAL.E5.OS.f_subcarrier = 15.345e6;

const.GAL.E5a.OS.f_carrier = 1176.45e6;
const.GAL.E5a.OS.T_code = 1e-3;
const.GAL.E5a.OS.L_code = 10230;
const.GAL.E5a.OS.f_chip = 15.345e6;
const.GAL.E5a.OS.f_subcarrier = 15.345e6;

const.GAL.E5b.OS.f_carrier = 1207.14e6;
const.GAL.E5b.OS.T_code = 1e-3;
const.GAL.E5b.OS.L_code = 10230;
const.GAL.E5b.OS.f_chip = 15.345e6;
const.GAL.E5b.OS.f_subcarrier = 15.345e6;

% Source: Hann01, p. 37
const.conductivity.concrete = 2e-5; % Sigma, conductivity
const.conductivity.dryGround = 1e-5;
const.conductivity.medDryGround = 4e-2;
const.conductivity.wetGround = 2e-1;
const.conductivity.freshWater = 2e-1;
const.conductivity.seaWater = 4;
const.conductivity.glass = 0.0025;
const.conductivity.steel = 1.45e6;

const.relativePermittivity.concrete = 3; % Epsilon, relative permittivity
const.relativePermittivity.dryGround = 4;
const.relativePermittivity.medDryGround = 7;
const.relativePermittivity.wetGround = 30;
const.relativePermittivity.freshWater = 80;
const.relativePermittivity.seaWater = 20;
const.relativePermittivity.glass = 7;
const.relativePermittivity.steel = 1;