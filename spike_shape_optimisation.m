%% Script to determine the optimal shape of an Aerospike nozzle
%Primary resource: (1) Rei Masuda. Optimal design of annular aerospike engine nozzle. Ann Arbor: California State University, Long Beach; 2002.

clc 
clear

load sim_results.mat
load universalConstants.mat
load rocketDesignParams.mat
load InitialConfigVars.mat
load rocketConfig.mat

n_points = 100;

tol = 0.0001;

[T_flame, gamma, m_mol, R,c_star] = thermochem(OF,P_cc,etac);

mdot_prop = m_prop_total/t_burn_result ;%average propellant mass flow

A_throat_check = (c_star*mdot_prop)/(P_cc*g0);

if abs(A_throat_check - A_throat) > tol
    
    disp('Error: Throat area mismatch. Consistency of analysis cannot be guaranteed')
end

A_throat = A_throat_check;

P_throat = P_cc*(1 + (gamma-1)/2)^(gamma/(gamma-1)); %subscript is unclear, check this

M_exit = sqrt((((P_throat/P_amb)^((gamma-1)/gamma))-1)*(2/(gamma-1))); %exit mach number of the flow

M = linspace(1,M_exit,n_points);

theta_throat = sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*(M_exit^2 - 1))) - atan(sqrt(M_exit^2 - 1)); %effectively nu(M_exit), the PM Angle for M_exit

A_exit = A_throat * (M_exit*((2/(gamma+1))*(1+((gamma-1)/2)*(M_exit^2)))^((-gamma+1)/(2*(gamma-1))))^(-1); %effective exit area of spike

r_exit = sqrt(A_exit/pi);

epsilon = A_exit/A_throat;

for i = 1 : n_points
    nu(i) = sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*(M(i)^2 - 1))) - atan(sqrt(M(i)^2 - 1)); %general PM angle equation

    mu(i) = asin(1/M(i)); %mach angle equation
    
    alpha(i) = theta_throat - nu(i) + mu(i); %arbitraty relation for clarity
    
    beta(i) = (1/M(i))*((2/gamma+1)*(1+(gamma-1/2)*(M(i)^2)))^((gamma-1)/2*(gamma-1)); %equal to A/At, inverse isentropic area ratio
    
    r(i) = r_exit*sqrt(1 - ((beta(i)/epsilon)*cos(theta_throat-nu(i))));

    l(i) = (r_exit-r(i))/sin(alpha(i));
end

l_tot = (r_exit-r(n_points)/sin(alpha(n_points)))*cos(alpha(n_points));

x = linspace(1,l_tot,n_points);

plot(x,r)