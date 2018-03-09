%% Script to determine the optimal shape of an Aerospike nozzle
%Primary resource: (1) Rei Masuda. Optimal design of annular aerospike engine nozzle. Ann Arbor: California State University, Long Beach; 2002.

clc 
clear

load sim_results.mat
load universalConstants.mat
load rocketDesignParams.mat

[T_flame, gamma, m_mol, R,c_star] = thermochem(OF,P_cc,etac);

A_thoat_check = (c_star*mdot_prop)/(P_cc*g);

if abs(A_throat_check - A_throat) > tol
    
    disp('Error: Throat area mismatch. Consistency of analysis cannot be guaranteed')
end

P_throat %subscript is unclear, check this

M_exit = sqrt((((P_throat/P_atmo)^((gamma-1)/gamma))-1)*(2/(gamma-1))); %exit mach number of the flow

nu = sqrt((gamma+1)/(gamma-1))*arctan(sqrt(((gamma-1)/(gamma+1))*(M^2 - 1))) - arctan(sqrt(M^2 - 1)); %general PM angle equation

mu = arcsin(1/M); %mach angle equation

theta_throat = sqrt((gamma+1)/(gamma-1))*arctan(sqrt(((gamma-1)/(gamma+1))*(M_exit^2 - 1))) - arctan(sqrt(M_exit^2 - 1)); %effectively nu(M_exit), the PM Angle for M_exit

alpha = theta_throat - nu + mu; %arbitraty relation for clarity

beta = (1/M)*((2/gamma+1)*(1+(gamma-1/2)*(M^2)))^((gamma-1)/2*(gamma-1)); %equal to A/At, inverse isentropic area ratio

A_exit = A_throat*(1/M_exit)*((2/gamma+1)*(1+(gamma-1/2)*(M_exit^2)))^((gamma-1)/2*(gamma-1)); %effective exit area of spike



