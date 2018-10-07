function rocketDesign = intialConfig(universalConstants, rocketDesignParameters, InitialConfigVars, regRateParams)


%% extract variables so we can use them easily

%Universal Constants
g0    = universalConstants.g0;
bar   = universalConstants.bar;
P_amb = universalConstants.P_amb;

%rocketDesignParameters

if rocketDesignParameters.port.type=="Dport"
    D_outer              = rocketDesignParameters.port.initialGuess.D_outer;
    tau                  = rocketDesignParameters.port.initialGuess.tau;
    fuelweb_initialguess = rocketDesignParameters.port.initialGuess.fuelweb_initialguess;
end


I_total= rocketDesignParameters.I_total;
F_init = rocketDesignParameters.F_init;
t_burn  = I_total/F_init;

P_cc     = rocketDesignParameters.P_cc;
lambda   = rocketDesignParameters.lambda;
rho_fuel = rocketDesignParameters.rho_fuel;
etac     = rocketDesignParameters.etac;
GO_max   = rocketDesignParameters.GO_max;

%InitialConfigVars
OF       = InitialConfigVars.OF;
T_req    = InitialConfigVars.T_req;
Isp_init = InitialConfigVars.Isp_init;
Isp_avg  = InitialConfigVars.Isp_avg;


%regRate
a = regRateParams.a;
n = regRateParams.n;
m = regRateParams.m;

%% Start configuration


%% Choose a OF, Then Mass and mass flows

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

[P_vap, rho_ox, dens_vap] = nitrous(T_req);

P_vap = P_vap*bar; %convert to Pascals

mdot_propinit = F_init/(Isp_init*g0);    %initial mass flow rate (SPAD eq 7.79)
mdot_fuelinit = mdot_propinit/(1+OF); %[SPAD, eq 7.79]
mdot_oxinit = mdot_propinit-mdot_fuelinit; %[SPAD, eq 7.79]

%% Determine C* Cstar [PROPEP?]


%This function calculates the values below.
[T_flame, gamma, m_mol, R,c_star] = thermochem(OF,P_cc,etac);

%gamma = 1.24; % ratio of specific heats (guess, but should be properly determined based on chosen O/F)
%m_mol = 0.0262109; %Molar mass (kg/mol)
%R = 8.314/m_mol; %Specific Gas Constant [SPAD, eq 7 .72]
%T_flame = 3300; %[K] Guess!! check with [Propep]
%c_star = etac*sqrt(gamma*R*T_flame)/(gamma*(2/(gamma+1))^((gamma+1)/(2*gamma-2))); %characteristic velocity [SPAD, eq 7.71]


%% Determine pressure levels

%Combustion Chamber pressure cannot be higher than Max tank pressure.


%dPvalve = ((mdoto/Cdis_valve*Avalve)^2)*1/(2*rhoo); %required pressure drop over valve [RPE Page 282]
%Pcc + dPinj + dPvalve = Ptank
%Note Cdis = mdot_actual/mdot_theoretical
Cdis_inj = 0.9; %[RPE page 280]
Cdis_valve = 1; %[guess]
K_inj = 1.5; %guess -> determine with water for preferred injector type

%P_vap = 55*bar; %[bar] room temp vapour pressure for nitrous, [physics of nitrous oxide] in the drive
P_tank = P_vap; %given, "nitrous" (?) [see engine_config_v6
%P_cc = Pcc_max(A_throat,mdot_prop,c_star); %max chamber pressure limited by tank pressure, materials, etc.


dP_inj = 0.15*P_cc; %must be high enough to isolate feed system from pressure transients in the CC [SPAD Page 232]
A_inj = mdot_propinit*sqrt(K_inj/(2*rho_fuel*dP_inj));
d_inj = 2*sqrt(A_inj/pi); %diamter of 1 injector hole
drill_lim_inj = 0.5e-3; %[m] limit of hole diameter for injector
holenum_inj = d_inj/drill_lim_inj;
%A_valve = 0.02; %standin valu

%A_inj = (mdot_oxinit^2/(2*rho_ox*Cdis_inj^2))*(P_tank - P_cc - (mdot_oxinit/(Cdis_valve*A_valve))^2*(1/(2*rho_ox)))^(-1); % WRONG FIX IT    Design output for area of injector orifices
%dP_inj = ((mdot_oxinit/Cdis_inj*A_inj)^2)*1/(2*rho_ox); %required pressure drop over injector.


%% configure combustion port
%assume single cylindrical port

GO_init = GO_max; %[kg/(m^2*s)] initial oxidiser flow flux [SPAD says blow-off/flooding limit is usually at 350-700 kg/m^2 s, so we used 350 to be safe]

A_port = mdot_oxinit/GO_init; %[m^2] [SPAD, eq. 7.82]

rocketDesign.port.type = rocketDesignParameters.port.type;


switch rocketDesignParameters.port.type
    case "circ"
        Diameter_port_init = 2*sqrt(A_port/pi); %[m] initial port diameter
        Perimeter_port = pi*Diameter_port_init;
        
        rocketDesign.port.IntialDiameter   = Diameter_port_init;
        rocketDesign.port.InitialPerimeter = Perimeter_port;
        
    case "Dport"
        %need to solve for fuel web, we have assumed values for D_outer and
        %for tau in dportInitialGuesses.mat
        
        fun = @(fuelweb) (DPort(D_outer,fuelweb,tau)-A_port);
        fuelweb_initial=fzero(fun, fuelweb_initialguess); %outputs the required fuelweb thickness
        %check if outputs make sense
        if 2*fuelweb_initial>D_outer-tau
            disp('Error in creating initial fuel grain: chamber area is negative')
        end
        %should perform a human check of whether the numbers are reasonable.
        PortParameters = [D_outer,fuelweb_initial,tau]; %vector of port parameters
        [~,Perimeter_port] = DPort(D_outer,fuelweb_initial,tau);
        
        rocketDesign.port.D_outer = D_outer;
        rocketDesign.port.InitialFuelWeb = fuelweb_initial;
        rocketDesign.port.tau = tau;
        rocketDesign.port.InitialPerimeter = Perimeter_port;
end
%Should perform a human check on the suitability of these numbers


GF_init = GO_init/OF; %initial fuel flow flux [SPAD, eq 7.83]

%regression rate formula takes form of:
% r= a G^n L^m where the coefficients are defined in regRateParams.mat

Lp = (mdot_fuelinit/(rho_fuel*a*(GO_init+GF_init)^n*Perimeter_port))^(1/(m+1)); %length of port (m) [SPAD, eq 7.88]


if rocketDesignParameters.port.type == "circ"
    %only for circular:
    Diameter_port_fin = sqrt((4*m_f)/(pi*Lp*rho_fuel)+Diameter_port_init^2);   %final diameter of the port [SPAD, eq 7.95]
    
    fuelweb=(Diameter_port_fin-Diameter_port_init)/2; %thickness of fuel that gets burnt [SPAD, 7.96]
    
    PortParameters = [Diameter_port_fin,fuelweb];
    
    rocketDesign.port.FinalDiameter=Diameter_port_fin;
    rocketDesign.port.InitialFuelWeb=fuelweb;
end

r = a*((GO_init+GF_init)^n)*(Lp^m); %calculate reg rate for reference



% N Port Wagon Wheel: ALL FROM SPAD FIG. 7.23
%n_max = 3;
%wagon_config = wagon_wheel_geometry(n_max,A_port,mdot_fuelinit,GO_init,a,m,n,rho_fuel,m_f);


%% determine nozzle area

%Since the combustion chamber diameter (dport final) and pressure is known,
%we can determine the throat area required to choke the flow.

A_throat = mdot_propinit*etac*sqrt(gamma*R*T_flame)*((gamma+1)/2)^((gamma+1)/(2*gamma-2))/(P_cc*gamma);


%note for above: used the T_flame as total temperature of flow which is
%probably inaccurate.

%%% Use F = mdot Ve = mdot Me sqrt(gamma R T0/(1+(gamma-1)/2*Me^2)) and
%%% then solve for Me

syms Me;
M_exit_target=double(vpasolve(mdot_propinit*Me*sqrt(gamma*R*T_flame/(1+(gamma-1)/2*Me^2)) == F_init, Me, 3)); %solve for required exit mach for desired thrust

[~, ~, ~, ~, expansionRatio] = flowisentropic(gamma, M_exit_target,'mach'); %this is a function in the matlab aerospace toolbox

A_exit = expansionRatio*A_throat;


%% key outputs:
rocketDesign.intialRegRate=r;

rocketDesign.c_star = c_star;

rocketDesign.m_f = m_f;

rocketDesign.m_ox = m_ox;

rocketDesign.mdot_fuelinit = mdot_fuelinit;

rocketDesign.mdot_oxinit = mdot_oxinit;

rocketDesign.Lp = Lp;

rocketDesign.A_inj = A_inj;

rocketDesign.A_throat = A_throat;

rocketDesign.A_exit = A_exit;

rocketDesign.expansionRatio=expansionRatio;


%other design inputs that need to be carried over to rocketDesign
rocketDesign.etac = rocketDesignParameters.etac;

rocketDesign.lambda = rocketDesignParameters.lambda;

rocketDesign.rho_fuel = rocketDesignParameters.rho_fuel;

rocketDesign.P_cc=rocketDesignParameters.P_cc;

%%
%useful line for debugging
%plotCrossSection(porttype,PortParameters);

%% References:

%%% [SPAD]: Space Propulsion Analysis and Design (Humble, book)
%%% [2]: DEVELOPMENT OF SCALABLE SPACE-TIME AVERAGED REGRESSION RATE
%%% EXPRESSIONS FOR HYBRID ROCKETS,  by M. Arif Karabeyoglu, Brian J.
%%% Cantwell and Greg Zilliac (AIAA 2005-3544)



end
