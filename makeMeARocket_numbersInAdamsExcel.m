clear all; %close all; clc;

%%% This script defines the parameters that need to be defined and then
%%% runs all necessary subscripts to be able to output and plot the
%%% performance of the rocket with time. This script should make it very
%%% easy to modify input parameters to observe their effects on the output
%%% parameters so we can tweak the performance to match that which is
%%% desired.

%%% First created by Dev on 23 Feb 2018

%% Inputs

%Universal Constants
g0      = 9.81;         %[m/s^2]
bar     = 100000;       %[Pa]
P_amb   = 1.0135*bar;   %[Pa] ambient pressure;

save universalConstants.mat g0 bar P_amb


%Rocket Design Parameters and targets
porttype = 1;
%1 is for cylinder
%2 is for D port

I_total = 2000;      %[Ns] Input goal total impulse, choose considering Adam Baker's Tank
F_init  = 500;       %[N] Input Goal average thrust
t_burn  = I_total/F_init; %[s] total burn time in seconds

P_cc     = 30*bar; %Chosen based on limits of tank pressurisation, and recommendations of [Physics of nitrous oxide, Aspire Space] in the drive
lambda   = 0.85;      %nozzle thrust efficiency factor
rho_fuel = 953;    %density of fuel (kg/m^3), guess, should be determined more accurately
etac     = 0.95;   %combustion efficiency [SPAD]
GO_max   = 350;    %[kg/(m^2*s)] max oxidiser flow flux [SPAD says blow-off/flooding limit is usually at 350-700 kg/m^2 s, so we used 350 to be safe]

deltaT=0.01; %[s] timestep in time simulator script
%performed a quick test to see that if deltaT is made smaller by one order
%of magnitude the results are unaffected (which is awesome - by dev, 23 Feb
%2018, 10:10pm)

save rocketDesignParams.mat porttype I_total F_init t_burn P_cc lambda rho_fuel etac GO_max deltaT


%Initial guesses (as defined for initialConfig)
OF       = 7;   %intial oxidiser to fuel ratio (guess used for config)
T_req    = 25;   %[degrees C] Input for 'nitrous;' function: temperature of ox tank contents
Isp_init = 200; %initial guess, should be iterated to maximise overall average Isp
Isp_avg  = 0.98*Isp_init;   %[SPAD 7.4.4] Says that the average Isp should be ~2.0% lower than initial


if porttype == 2
    D_outer = 80e-3; %[m]
    tau     = 2e-3;  %[m] central metal half-thickness
    fuelweb_initialguess = 10e-3; %[m] initial guess -> use this to control output
    save dportInitialGuesses.mat D_outer tau fuelweb_initialguess
end

save InitialConfigVars.mat  OF Isp_init Isp_avg T_req


%Regression Rate Params
%regression rate formula takes form of:
% r= a G^n L^m where
a = 2.34e-5;  %need correct values! Using Paraffin/GOX value from [ref 2, Table 1] - gives rdot in m/sec but uses SI units for
n = 0.8;      %need correct values!
m = -0.2;     %need correct values!

%using numbers from adam bakers excel file.
a = 0.0000236;
n = 0.605;
m = 0;




save regRateParams.mat a n m

%outputConditions
qplot = 1; %if 1, there will be an output displayed
qdisp = 1; %if 1, there will be an output displayed
save outputConditions.mat qplot qdisp

%notes: Determine pressure levels variables have not been transported over

%%
%******************
% NOTHING NEEDS TO BE CHANGED BELOW THIS FOR NORMAL USAGE
%******************


%% Run Config script

engineConfig_Initial_2018;

%% Peform Time based simulation

timeSimulator

%% Output results

load outputConditions.mat

if porttype==1
    Dia_store=PortParameters_store(:,1)-2*PortParameters_store(:,2);
end

clc;
disp('Simulation Successful')

if qdisp == 1
    disp('All values are in SI units')
    I_total_result
    t_burn_result
    F_init_result
    F_avg_result
    Isp_avg
    m_f_total
    m_ox_total
    m_prop_total
end
if qplot == 1
    
    t_axis = 0:deltaT:t_burn_result;
    
    figure(1);
    
    
    subplot(3,4,1)
    plot(t_axis,F,[0],[0])
    title('F')
    
    subplot(3,4,2)
    plot(t_axis,c_star,[0],[0])
    title('c star')
    
    
    subplot(3,4,3)
    plot(t_axis,PortParameters_store(:, 2),[0],[0])
    
    title('Fuel web thickness [m]')
    
    subplot(3,4,4)
    plot(t_axis,G_ox,[0],[0])
    title('G ox')
    
    subplot(3,4,5)
    plot(t_axis,G_fuel,[0],[0])
    title('G fuel')
    
    subplot(3,4,6)
    plot(t_axis,G_prop,[0],[0])
    title('G prop')
    
    subplot(3,4,7)
    plot(t_axis,mdot_fuel,[0],[0])
    title('mdot fuel')
    
    subplot(3,4,8)
    plot(t_axis,Isp,[0],[0])
    title('Isp')
    
    subplot(3,4,9)
    plot(t_axis,M_exit,[0],[0])
    title('M exit')
    
    subplot(3,4,10)
    plot(t_axis,OF,[0],[0])
    title('OF')
    
    subplot(3,4,11)
    plot(t_axis,P_cc,[0],[0])
    title('P cc')
    
    subplot(3,4,12)
    plot(t_axis,P_exit,[0],[0])
    title('P exit')
    
    % Plot cross sections, before and after
    figure (2);
    clf
    
    subplot(1, 2, 1);
    plotCrossSection(porttype, PortParameters_store(1, :));
    title('Intial cross section [m]')
    
    subplot(1, 2, 2);
    plotCrossSection(porttype, PortParameters_store(end, :));
    title('Final cross section [m]')
end

