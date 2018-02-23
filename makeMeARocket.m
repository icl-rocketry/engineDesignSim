clear all; close all; clc;

%%% This script defines the parameters that need to be defined and then
%%% runs all necessary subscripts to be able to output and plot the
%%% performance of the rocket with time. This script should make it very
%%% easy to modify input parameters to observe their effects on the output
%%% parameters so we can tweak the performance to match that which is
%%% desired.

%%% First created by Dev on 23 Feb 2018

%% Inputs

%Universal Constants
g0 = 9.81; %[m/s^2]
bar = 100000; %[Pa]
P_amb = 1.0135*bar;  %[Pa] ambient pressure;

save universalConstants.mat g0 bar P_amb

%Rocket Design Parameters and targets
porttype = 2;

I_total = 0.55*6000; %[Ns] Input goal total impulse, choose considering Adam Baker's Tank
F_init = 0.85*600; %[N] Input Goal average thrust
t_burn = I_total/F_init; %[s] total burn time in seconds


P_cc = 25*bar;  %Chosen based on limits of tank pressurisation, and recommendations of [Physics of nitrous oxide, Aspire Space] in the drive
lambda = 1;     %nozzle thrust efficiency factor
rho_fuel=953; %density of fuel (kg/m^3), guess, should be determined more accurately
etac = 0.95; %combustion efficiency [SPAD]
GO_max = 350; %[kg/(m^2*s)] max oxidiser flow flux [SPAD says blow-off/flooding limit is usually at 350-700 kg/m^2 s, so we used 350 to be safe]


save rocketDesignParams.mat porttype I_total F_init t_burn P_cc lambda rho_fuel etac GO_max







%Initial guesses (as defined for initialConfig)

OF = 3; %intial oxidiser to fuel ratio (guess used for config)
Isp_init=200; %initial guess, should be iterated to maximise overall average Isp
Isp_avg = 0.98*Isp_init; %[SPAD 7.4.4] Says that the average Isp should be ~2.0% lower than initial
T_req = 5; %[degrees C] Input for 'nitrous;' function: temperature of ox tank contents

if porttype == 2
    D_outer = 80e-3; %[m]
    tau = 2e-3; %[m] central metal half-thickness
    fuelweb_initialguess = 10e-3;%[m] initial guess -> use this to control output
    save dportInitialGuesses.mat D_outer tau fuelweb_initialguess
end

save InitialConfigVars.mat  OF Isp_init Isp_avg T_req 


%Regression Rate Params
%regression rate formula takes form of:
% r= a G^n L^m where
a=2.34e-5; %need correct values! Using Paraffin/GOX value from [ref 2, Table 1] - gives rdot in m/sec but uses SI units for
n=0.8; %need correct values!
m=-0.2; %need correct values!

save regRateParams.mat a n m


%notes: Determine pressure levels variables have not been transported over

%% Run Config script

engineConfig_Initial_2018;

%% Peform Time based simulation

%% Output results