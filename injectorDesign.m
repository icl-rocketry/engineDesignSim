%% Injector Design
% Will Harradence
% Imperial Aeronautics 2019/20

% Source is "Modern Engineering for the Design of Liquid Propellant Rocket
% Engines", by Huzel and Huang, P.113

%% Setup
clear
clc

load n2o.mat
bar = 1e5;

Pcc = 35*bar;
dP = 0.2*Pcc; % 20% chamber pressure drop across injector
T = 300; %K
n2o.T = T;
rho = densityLiquid(n2o);

k = 1.7; % Head loss factor, sharp-edged orifices
%k = 1.3; %radiused orifices

m_dot = [0.1 0.2 0.3 0.4]; %Mass flow rates desired, kg/s
n = 3; %number of orifices desired

%% Unit Conversion
m_doti = convmass(m_dot,'kg','lbm');
dPi = convpres(dP,'Pa','psi');
rhoi = convdensity(rho,'lbm/ft^3','kg/m^3');

%% Area Calculation
%Gives total injector area required

Ai = m_doti.*sqrt((2.238*k)/(rhoi*dPi)); %result in in^2

A = 0.00064516*Ai; % convert from in^2 to m^2
A_small = A.*1000000; %mm^2

%% Diameter Calculation

di = ((3.627*k.*m_doti.^2)./(rhoi*dPi.*n.^2)).^0.25; %result in inches

d = convlength(di,'in','m'); % convert from in to m
d_small = d.*1000; %mm

%% Injector Velocity Calculation
%WARNING: little information is given on this calculation in the source,
%proceed with caution

Vi = ((144.*m_doti)./(Ai.*rhoi)); % ft/s

V = convvel(Vi,'ft/s','m/s');

