function [injectorDesign] = injectorConfig(injectorDesignParams,rocketDesign)
%% Injector Design
% Will Harradence
% Imperial Aeronautics 2019/20

% Source is "Modern Engineering for the Design of Liquid Propellant Rocket
% Engines", by Huzel and Huang, P.113

%% Unpacking

n = injectorDesignParams.n;
throttle_states = injectorDesignParams.throttle_states;
m_dot_nominal = rocketDesign.mdot_oxinit;
orifice = injectorDesignParams.orificeType;
P_cc = rocketDesign.P_cc;

dP = 0.2*P_cc; % 20% chamber pressure drop across injector

rho = 725.3193;

%% Setting head loss factor
if orifice == 'sharp'
    
    k = 1.7; % Head loss factor, sharp-edged orifices

elseif orifice == 'radiused'
    
    k = 1.3;
    
else
    
    disp('Improper injector orifice type selection, please check spelling and supported orifice types')
    return
    
end

%% Throttling limits
% Sets the mass flow rates for different specified throttle states.
for i  = 1 : length(throttle_states)
    
    m_dot(i) = m_dot_nominal*throttle_states(i);

end

%% Unit Conversion
m_doti = convmass(m_dot,'kg','lbm');
dPi = convpres(dP,'Pa','psi');
rhoi = convdensity(rho,'lbm/ft^3','kg/m^3');

%% Area Calculation
% Gives total injector area required, across all orifices

Ai = m_doti.*sqrt((2.238*k)/(rhoi*dPi)); %result in in^2

A = 0.00064516*Ai; % convert from in^2 to m^2
A_small = A.*1000000; %mm^2 , good for debug

%% Diameter Calculation
% Gives the orifice diameter required to achieve the correct pressure drop
% across the injector.

di = ((3.627*k.*m_doti.^2)./(rhoi*dPi.*n.^2)).^0.25; %result in inches

d = convlength(di,'in','m'); % convert from in to m
d_small = d.*1000; %mm, good for debug

%% Injector Velocity Calculation
% WARNING: little information is given on this calculation in the source,
% proceed with caution

Vi = ((144.*m_doti)./(Ai.*rhoi)); % ft/s

V = convvel(Vi,'ft/s','m/s');

%% Set Outputs
injectorDesign.A_inj  = A;
injectorDesign.d_inj = d;
injectorDesign.V_inj = V;

end

