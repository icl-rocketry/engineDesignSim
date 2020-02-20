function rocketDesign = intialConfig(universalConstants, rocketDesignParameters, InitialConfigVars, regRateParams)
%% Unpacking

%Universal Constants
g0    = universalConstants.g0;
bar   = universalConstants.bar;
P_amb = universalConstants.P_amb;

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


%% Calculate Flow rates
[P_vap, rho_ox, dens_vap] = nitrous(T_req);

P_vap = P_vap*bar; %convert to Pascals

mdot_propinit = F_init/(Isp_init*g0);    %initial mass flow rate (SPAD eq 7.79)
mdot_fuelinit = mdot_propinit/(1+OF); %[SPAD, eq 7.79]
mdot_oxinit = mdot_propinit-mdot_fuelinit; %[SPAD, eq 7.79]

%% Determine C* Cstar [PROPEP?]


%This function calculates the values below.
%[T_flame, gamma, m_mol, R,c_star] = thermochem(OF,P_cc,etac);
[T_flame, gamma, m_mol, R,c_star] = paraffinThermo(etac);

%gamma = 1.24; % ratio of specific heats (guess, but should be properly determined based on chosen O/F)
%m_mol = 0.0262109; %Molar mass (kg/mol)
%R = 8.314/m_mol; %Specific Gas Constant [SPAD, eq 7 .72]
%T_flame = 3300; %[K] Guess!! check with [Propep]
%c_star = etac*sqrt(gamma*R*T_flame)/(gamma*(2/(gamma+1))^((gamma+1)/(2*gamma-2))); %characteristic velocity [SPAD, eq 7.71]


%% Configure Combustion Port
%assume single cylindrical port

GO_init = GO_max; %[kg/(m^2*s)] initial oxidiser flow flux [SPAD says blow-off/flooding limit is usually at 350-700 kg/m^2 s, so we used 350 to be safe]

A_port = mdot_oxinit/GO_init; %[m^2] [SPAD, eq. 7.82]

Diameter_port_init = 2*sqrt(A_port/pi); %[m] initial port diameter
Perimeter_port = pi*Diameter_port_init;

rocketDesign.port.IntialDiameter   = Diameter_port_init;
rocketDesign.port.InitialPerimeter = Perimeter_port;
        
%Should perform a human check on the suitability of these numbers

GF_init = GO_init/OF; %initial fuel flow flux [SPAD, eq 7.83]

%regression rate formula takes form of:
% r= a G^n L^m where the coefficients are defined in regRateParams.mat

Lp = (mdot_fuelinit/(rho_fuel*a*(GO_init+GF_init)^n*Perimeter_port))^(1/(m+1)); %length of port (m) [SPAD, eq 7.88]


%only for circular:

Diameter_port_fin = sqrt((4*m_f)/(pi*Lp*rho_fuel)+Diameter_port_init^2);   %final diameter of the port [SPAD, eq 7.95]

fuelweb=(Diameter_port_fin-Diameter_port_init)/2; %thickness of fuel that gets burnt [SPAD, 7.96]

PortParameters = [Diameter_port_fin,fuelweb];

r = a*((GO_init+GF_init)^n)*(Lp^m); %calculate reg rate for reference

%% determine nozzle area

%Since the combustion chamber diameter (dport final) and pressure is known,
%we can determine the throat area required to choke the flow.

A_throat = mdot_propinit*etac*sqrt(gamma*R*T_flame)*((gamma+1)/2)^((gamma+1)/(2*gamma-2))/(P_cc*gamma);

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

rocketDesign.A_throat = A_throat;

rocketDesign.A_exit = A_exit;

rocketDesign.expansionRatio=expansionRatio;

rocketDesign.port.FinalDiameter=Diameter_port_fin;

rocketDesign.port.InitialFuelWeb=fuelweb;

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
