clear all; 
clc;
close all;

%% Version 1
% Last edited by Dev & Will on Feb 9 2018 1801

%%% This script aims to configure the rocket using Chapter 7 of Space
%%% Propulsion Analysis and Design (SPAD)
%%% Rocket Propulsion Elements (RPE)

%% Constants
g0 = 9.81; %[m/s^2]
bar = 100000; %[Pa] 
%% Design Inputs/Constants

I_total = 6000; %[Ns] Input goal total impulse, choose considering Adam Baker's Tank
F_init = 600; %[N] Input Goal average thrust
t_burn = I_total/F_init; %[s] total burn time in seconds


%% Choose a OF, Then Mass and mass flows 

OF = 7; %oxidiser to fuel ratio (guess, but should be determined to maximise average Isp)

Isp_init=200; %initial guess, should be iteratively maximised
Isp_avg = 0.98*Isp_init; %[SPAD 7.4.4] Says that the average Isp should be ~2.0%lower than initial
% Isp = I0/mprop*g0
m_prop = I_total/(Isp_avg*g0); %[kg] RPE Ch.2
m_ox = m_prop*OF/(OF+1);
m_f = m_prop - m_ox;


%% size fuel tanks

%there are two methods: either 
% (A) sizing using mission anlysis:
%   (1) mpayload
%   (2) Isp:    take avgIsp=0.99*maxIsp
%   (3) delv:   mission analysis tells you this
%   (4) finert: around 0.16-0.20 for solid, a bit lower for hybirds
% or 
% (B) sizing for target performance:
%   (1) It target (total impulse target)
%   (2) Isp
%   (3) target thrust

T_req = 5; %[degrees C] Input for 'nitrous;' function: temperature of ox tank contents
[P_vap, rho_ox, dens_vap] = nitrous(T_req);
rho_f=953; %density of fuel (kg/m^3), guess, should be determined more accurately
P_vap = P_vap*bar;

mdot_propinit = F_init/(Isp_init*g0);    %initial mass flow rate (SPAD eq 7.79)
mdot_fuelinit = mdot_propinit/(1+OF); %[SPAD, eq 7.79]
mdot_oxinit = mdot_propinit-mdot_fuelinit; %[SPAD, eq 7.79]


%% Determine C* Cstar [PROPEP?]

etac = 0.95; %combustion efficiency [SPAD]
gamma = 1.24; % ratio of specific heats (guess, but should be properly determined based on chosen O/F)
m_mol = 0.0262109; %Molar mass (kg/mol)
R = 8.314/m_mol; %Specific Gas Constant [SPAD, eq 7 .72]
T_flame = 3300; %[K] Guess!! check with [Propep]
c_star = etac*sqrt(gamma*R*T_flame)/(gamma*(2/(gamma+1))^((gamma+1)/(2*gamma-2))); %characteristic velocity [SPAD, eq 7.71]


%% Determine pressure levels

%Combustion Chamber pressure cannot be higher than Max tank pressure. 


%dPvalve = ((mdoto/Cdis_valve*Avalve)^2)*1/(2*rhoo); %required pressure drop over valve [RPE Page 282]
%Pcc + dPinj + dPvalve = Ptank
%Note Cdis = mdot_actual/mdot_theoretical
Cdis_inj = 0.9; %[RPE page 280]
Cdis_valve = 1; %[guess]

%P_vap = 55*bar; %[bar] room temp vapour pressure for nitrous, [physics of nitrous oxide] in the drive
P_tank = P_vap; %given, "nitrous" (?) [see engine_config_v6
%P_cc = Pcc_max(A_throat,mdot_prop,c_star); %max chamber pressure limited by tank pressure, materials, etc.
P_cc = 25*bar; %Chosen based on limits of tank pressurisation, and recommendations of [Physics of nitrous oxide, Aspire Space] in the drive

A_valve = 0.02; %standin value

A_inj = (mdot_oxinit^2/(2*rho_ox*Cdis_inj^2))*(P_tank - P_cc - (mdot_oxinit/(Cdis_valve*A_valve))^2*(1/(2*rho_ox)))^(-1); % WRONG FIX IT    Design output for area of injector orifices
dP_inj = ((mdot_oxinit/Cdis_inj*A_inj)^2)*1/(2*rho_ox); %required pressure drop over injector. 

%% configure combustion port
%assume single cylindrical port

GO_init = 350; %[kg/(m^2*s)]initial oxidiser flow flux [SPAD says blow-off/flooding limit is usually at 350-700 kg/m^2 s, so we used 350 to be safe]

A_port = mdot_oxinit/GO_init; %[m^2] [SPAD, eq. 7.82]

Dport_init = 2*sqrt(A_port/pi); %[m] initial port diameter

GF_init = GO_init/OF; %initial fuel flow flux [SPAD, eq 7.83]

%regression rate formula takes form of:
% r= a G^n L^m where
a=2.34e-5; %need correct values! Using Paraffin/GOX value from [ref 2, Table 1] - gives rdot in m/sec but uses SI units for 
n=0.8; %need correct values!
m=-0.2; %need correct values!

Lp = (((mdot_fuelinit*(pi^(n-1)))/(a*((4*mdot_propinit)^n)*rho_f))*((Dport_init)^(2*n-1)))^(1/(m+1));%length of port (m) [SPAD, eq 7.91]

Dport_fin = sqrt((4*m_f)/(pi*Lp*rho_f)+Dport_init^2);   %final diameter of the port [SPAD, eq 7.95]

fuelweb=(Dport_fin-Dport_init)/2; %thickness of fuel that gets burnt [SPAD, 7.96]

r = a*((GO_init+GF_init)^n)*(Lp^m);
%% key outputs:
r

c_star

m_f 

m_ox

mdot_fuelinit

mdot_oxinit

Dport_init

Lp

fuelweb

A_inj

%% References:

%%% [SPAD]: Space Propulsion Analysis and Design (Humble, book)
%%% [2]: DEVELOPMENT OF SCALABLE SPACE-TIME AVERAGED REGRESSION RATE
%%% EXPRESSIONS FOR HYBRID ROCKETS,  by M. Arif Karabeyoglu, Brian J.
%%% Cantwell and Greg Zilliac (AIAA 2005-3544)

