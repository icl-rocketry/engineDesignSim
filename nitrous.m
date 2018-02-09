function [P, dens_liq, dens_vap] = nitrous(T_req)

%% VERSION 1
% Last edit by Joseph & Teejay
% 21:50 24/11/15

% Script based on p4 of The Physics of Nitrous Oxide (PNO) by Rick Newlands
% Input: Temperature T_req (degC)
% Output: Vapour pressure (P), liquid density (dens_liq), vapour density (dens_vap)

%% Data from PNO
T_data = [-20; -15; -10; -5; 0; 5; 10; 15; 20; 25; 30; 35]; %Temperature (degC)
P_data = [18.01; 20.83; 23.97; 27.44; 31.27; 35.47; 40.07; 45.1; 50.6; 56.6; 63.15; 70.33]; %Vapour pressure (bar)
dens_liq_data = [995.4; 975.2; 953.9; 931.4; 907.4; 881.6; 853.5; 822.2; 786.6; 743.9; 688; 589.4]; %Liquid density (kg/m^3)
dens_vap_data = [46.82; 54.47; 63.21; 73.26; 84.86; 98.41; 114.5; 133.9; 158.1; 190; 236.7; 330.4]; %Vapour density (kg/m^3)


%% Interpolation
[val i]  = min(abs(T_data-T_req)); %Index of closest value
stepT    = T_data(i+1)-T_data(i);  %Difference between temperatures in T_data
deltaT   = T_req-T_data(i);        %Difference between closest T and required T

%Outputs
P        = (P_data(i+1)-P_data(i))/stepT * deltaT + P_data(i); %Vapour pressure (bar)
dens_liq = (dens_liq_data(i+1)-dens_liq_data(i))/stepT * deltaT + dens_liq_data(i); %Liquid density (kg/m^3)
dens_vap = (dens_vap_data(i+1)-dens_vap_data(i))/stepT * deltaT + dens_vap_data(i); %Vapour density (kg/m^3)

end