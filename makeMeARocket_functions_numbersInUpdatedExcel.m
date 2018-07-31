%% make me a rocket with functions and structures

clear all; close all; clc;

%%

%Universal Constants

g0      = 9.81;         %[m/s^2]
bar     = 100000;       %[Pa]
P_amb   = 1.0135*bar;   %[Pa] ambient pressure;

universalConstants.g0      = g0;
universalConstants.bar     = bar;
universalConstants.P_amb   = P_amb;

%rocketDesignParameters

rocketDesignParameters.port.type='circ';

if rocketDesignParameters.port.type=="Dport"
    rocketDesignParameters.port.initialGuess.D_outer = 80e-3; %[m]
    rocketDesignParameters.port.initialGuess.tau     = 2e-3;  %[m] central metal half-thickness
    rocketDesignParameters.port.initialGuess.fuelweb_initialguess = 10e-3; %[m] initial guess -> use this to control output
end


rocketDesignParameters.I_total=1200;
rocketDesignParameters.F_init=180;
rocketDesignParameters.t_burn  = rocketDesignParameters.I_total/rocketDesignParameters.F_init; %[s] total burn time in seconds

rocketDesignParameters.P_cc     = 30*bar; %Chosen based on limits of tank pressurisation, and recommendations of [Physics of nitrous oxide, Aspire Space] in the drive
rocketDesignParameters.lambda   = 0.85;      %nozzle thrust efficiency factor
rocketDesignParameters.rho_fuel = 953;    %density of fuel (kg/m^3), guess, should be determined more accurately
rocketDesignParameters.etac     = 0.95;   %combustion efficiency [SPAD]
rocketDesignParameters.GO_max   = 350;    %[kg/(m^2*s)] max oxidiser flow flux [SPAD says blow-off/flooding limit is usually at 350-700 kg/m^2 s, so we used 350 to be safe]


%InitialConfigVars

InitialConfigVars.OF       = 7;   %intial oxidiser to fuel ratio (guess used for config)
InitialConfigVars.T_req    = 25;   %[degrees C] Input for 'nitrous;' function: temperature of ox tank contents
InitialConfigVars.Isp_init = 200; %initial guess, should be iterated to maximise overall average Isp
InitialConfigVars.Isp_avg  = 0.98*InitialConfigVars.Isp_init;   %[SPAD 7.4.4] Says that the average Isp should be ~2.0% lower than initial

%regRate
%regRateParams.a = 2.34e-5;  %need correct values! Using Paraffin/GOX value from [ref 2, Table 1] - gives rdot in m/sec but uses SI units for
%regRateParams.n = 0.8;      %need correct values!
%regRateParams.m = -0.2;     %need correct values!

%using numbers from adam bakers excel file.
regRateParams.a = 0.0000236;
regRateParams.n = 0.605;
regRateParams.m = 0;

%outputConditions
qplot = 1; %if 1, there will be an output displayed
qdisp = 1; %if 1, there will be an output displayed


%% run initial config

rocketDesign = intialConfig(universalConstants, rocketDesignParameters, InitialConfigVars, regRateParams);

%% multiply fuel web by k times

fuelwebFactorK=3;

%update fuel web
rocketDesign.port.InitialFuelWeb = rocketDesign.port.InitialFuelWeb*fuelwebFactorK;

%update outer diameter
switch rocketDesign.port.type
    case "circ"
        rocketDesign.port.FinalDiameter = rocketDesign.port.IntialDiameter+2*rocketDesign.port.InitialFuelWeb;
        
    case "Dport"
        error('NEED TO FIX THIS PART OF THE CODE')
        
end


%% run time based simulation

thresh=1e-6;
deltaT = 0.01;

[rocketPerformance, rocketSim] = timeBasedSimulator(universalConstants,rocketDesign,regRateParams, deltaT,thresh);


%% display results

disp('Simulation Successful')

if qdisp == 1
    disp('All values are in SI units')
    rocketDesign
    rocketPerformance
    avgRegRate_in_mm_per_s=1000*(rocketSim.port.fuelweb(1)-rocketSim.port.fuelweb(end))/rocketPerformance.t_burn_result
end


if qplot == 1
    
    t_axis = 0:deltaT:rocketPerformance.t_burn_result;
    
    figure(1);
    
    
    subplot(3,4,1)
    plot(t_axis,rocketSim.F,[0],[0])
    title('F')
    
    subplot(3,4,2)
    plot(t_axis,rocketSim.c_star,[0],[0])
    title('c star')
    
    
    subplot(3,4,3)
    plot(t_axis,rocketSim.port.fuelweb(1:end-1),[0],[0])
    
    title('Fuel web thickness [m]')
    
    subplot(3,4,4)
    plot(t_axis,rocketSim.G_ox,[0],[0])
    title('G ox')
    
    subplot(3,4,5)
    plot(t_axis,rocketSim.G_fuel,[0],[0])
    title('G fuel')
    
    subplot(3,4,6)
    plot(t_axis,rocketSim.G_prop,[0],[0])
    title('G prop')
    
    subplot(3,4,7)
    plot(t_axis,rocketSim.mdot_fuel,[0],[0])
    title('mdot fuel')
    
    subplot(3,4,8)
    plot(t_axis,rocketSim.Isp,[0],[0])
    title('Isp')
    
    subplot(3,4,9)
    plot(t_axis,rocketSim.M_exit,[0],[0])
    title('M exit')
    
    subplot(3,4,10)
    plot(t_axis,rocketSim.OF,[0],[0])
    title('OF')
    
    subplot(3,4,11)
    plot(t_axis,rocketSim.P_cc,[0],[0])
    title('P cc')
    
    subplot(3,4,12)
    plot(t_axis,rocketSim.P_exit,[0],[0])
    title('P exit')
    
    %Plot cross sections, before and after
    figure (2);
    clf
    
    
    subplot(1, 2, 1);
    
    plotPortXSection(rocketSim.port,1)
    
    title('Intial cross section [m]')
    
    subplot(1, 2, 2);
    plotPortXSection(rocketSim.port,length(t_axis))
    title('Final cross section [m]')
end









